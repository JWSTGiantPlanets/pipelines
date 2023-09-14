__version__ = '1.0.2'

from typing import Callable, TypeAlias

import numpy as np
import planetmapper
from astropy.io import fits
from scipy.ndimage import gaussian_filter

import tools

CoordinateTransform: TypeAlias = Callable[[int], int]
"""Transform from original to scaled pixel coordinates"""
ForwardModelFunc: TypeAlias = Callable[[planetmapper.BodyXY], np.ndarray]
"""Forward model function"""
FwhmFunc: TypeAlias = Callable[[float], float]
"""PSF FWHM function"""


# Note: for MIRI, PSF FWHM is roughly 3px for the longest wavelength of each channel


class PsfCorrectionError(Exception):
    """Base class for PSF correction errors."""


class TargetTooSmallError(PsfCorrectionError):
    """Raised when the target is too small to be scaled."""


class TargetNotInFovError(PsfCorrectionError):
    """Raised when none of the target is in the field of view."""


def correct_file(
    path: str,
    path_out: str | None = None,
    forward_model_func: ForwardModelFunc | None = None,
    psf_fwhm_func: FwhmFunc | None = None,
    mask_off_disc_pixels: bool = True,
    check_target_in_fov: bool = True,
    **kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Correct observed flux using a PSF model.

    The forward model (as produced by `forward_model_func`) is convolved with a Gaussian
    PSF (as defined by `psf_fwhm_func`) at each wavelength. The observed flux is then
    corrected by dividing by the PSF model. The 'SCI' and 'ERR' extensions of the
    observation are corrected.

    Args:
        path: Path to observation to correct.
        path_out: Path to save corrected observation to. If `None`, the corrected
            observation will not be saved.
        forward_model_func: Function to use to generate the forward model. This should
            take a `planetmapper.BodyXY` object as input and return a numpy array
            containing the forward model image. If `None`, the default function
            `default_forward_model_func` will be used. This function should generally
            have a value of 1 for illuminated pixels and 0 for unilluminated pixels.
        psf_fwhm_func: Function to use to calculate the PSF FWHM at each wavelength.
            This should take a wavelength in microns as input and return the PSF FWHM in
            arcseconds. If `None`, the default function `default_psf_fwhm_func` will be
            used.
        mask_off_disc_pixels: If `True`, pixels that are zero in the forward model will
            be set to NaN in the corrected observation.
        check_target_in_fov: If `True`, raise a `TargetNotInFovError` if none of the
            target is in the field of view.
        **kwargs: Additional keyword arguments to pass to
            `create_dynamically_scaled_body`.

    Returns:
        `(forward_model, psf_model_cube, cube_corrected)` tuple where `forward_model`
        is the forward model image, `psf_model_cube` is the PSF model cube, and
        `cube_corrected` is the corrected observation.

    Raises:
        TargetNotInFovError: If none of the target is in the field of view and
            `check_target_in_fov` is `True`.
    """
    if forward_model_func is None:
        forward_model_func = default_forward_model_func
    if psf_fwhm_func is None:
        psf_fwhm_func = default_psf_fwhm_func

    observation = planetmapper.Observation(path)

    forward_model = forward_model_func(observation)
    if check_target_in_fov and np.all(forward_model == 0):
        raise TargetNotInFovError('None of the target in field of view')

    body, transform, parameters = create_dynamically_scaled_body(observation, **kwargs)
    check_scaling(observation, body, transform)

    with fits.open(path) as hdul:
        data_header: fits.Header = hdul['SCI'].header  # type: ignore
        sci_cube: np.ndarray = hdul['SCI'].data  # type: ignore
        err_cube: np.ndarray = hdul['ERR'].data  # type: ignore
        wavelengths = tools.get_wavelengths(data_header)

        psf_model_cube = calculate_psf_modelled_cube(
            sci_cube.shape,
            wavelengths,
            body,
            transform,
            forward_model_func,
            psf_fwhm_func,
        )

        with tools.ignore_warnings('divide by zero encountered'):
            sci_cube_corrected = sci_cube / psf_model_cube
            err_cube_corrected = err_cube / psf_model_cube

        psf_fwhm_arcsec = np.array([psf_fwhm_func(wl) for wl in wavelengths])

        if mask_off_disc_pixels:
            sci_cube_corrected[:, forward_model == 0] = np.nan
            err_cube_corrected[:, forward_model == 0] = np.nan

        if path_out is not None:
            hdul['SCI'].data = sci_cube_corrected  # type: ignore
            hdul['ERR'].data = err_cube_corrected  # type: ignore

            header = fits.Header()
            header.add_comment('Forward model used for PSF correction.')
            header.add_comment(
                'This model is convolved with the PSF at each wavelength.'
            )
            header.add_comment('This is an unscaled version of the forward model.')
            hdul.append(
                fits.ImageHDU(
                    data=forward_model, header=header, name='PSF_FORWARD_MODEL'
                )
            )

            header = fits.Header()
            header.add_comment(
                'PSF Gaussian FWHM in arcseconds used for PSF correction.'
            )
            hdul.append(
                fits.ImageHDU(
                    data=psf_fwhm_arcsec, header=header, name='PSF_FWHM_ARCSEC'
                )
            )

            header = fits.Header()
            header.add_comment('PSF model cube used for PSF correction.')
            header.add_comment(
                'This is the scaled forward model convolved with the PSF at each wavelength.'
            )
            hdul.append(
                fits.ImageHDU(data=psf_model_cube, header=header, name='PSF_MODEL_CUBE')
            )

            header = hdul['PRIMARY'].header  # type: ignore
            tools.add_header_reduction_note(hdul, 'PSF corrected')
            header['HIERARCH PSF VERSION'] = (__version__, 'Software version')
            for k, v in parameters.items():
                k = f'HIERARCH PSF {k.upper()}'
                header[k] = v

            header['HIERARCH PSF FORWARD_MODEL'] = forward_model_func.__name__
            header['HIERARCH PSF PSF_FWHM'] = psf_fwhm_func.__name__
            header['HIERARCH PSF MASK_OFF_DISC_PIXELS'] = mask_off_disc_pixels
            header['HIERARCH PSF CHECK_TARGET_IN_FOV'] = check_target_in_fov

            tools.check_path(path_out)
            hdul.writeto(path_out, overwrite=True)

    return forward_model, psf_model_cube, sci_cube_corrected


def calculate_psf_modelled_cube(
    cube_shape: tuple[int, ...],
    wavelengths: np.ndarray,
    body_scaled: planetmapper.BodyXY,
    transform: CoordinateTransform,
    forward_model_func: ForwardModelFunc,
    psf_fwhm_func: FwhmFunc,
) -> np.ndarray:
    """
    Calculate a PSF modelled cube by convolving the forward model with a Gaussian PSF
    at each wavelength.

    Args:
        cube_shape: Shape of the cube to return.
        wavelengths: Wavelengths to calculate the PSF modelled cube at.
        body_scaled: Scaled version of the observation to forward model.
        transform: Transform from original to scaled pixel coordinates.
        forward_model_func: Function to use to generate the forward model. This should
            take a `planetmapper.BodyXY` object as input and return a numpy array
            containing the forward model image.
        psf_fwhm_func: Function to use to calculate the PSF FWHM at each wavelength.
            This should take a wavelength in microns as input and return the PSF FWHM in
            arcseconds.

    Returns:
        Numpy array containing the PSF modelled cube.
    """
    forward_model_scaled = forward_model_func(body_scaled)
    arcsec_per_px_scaled = body_scaled.get_plate_scale_arcsec()

    idxs1 = [transform(i) for i in range(cube_shape[1])]
    idxs2 = [transform(i) for i in range(cube_shape[2])]
    idxs = np.ix_(idxs1, idxs2)

    psf_cube = np.full(cube_shape, np.nan)
    for wl_idx, wl in enumerate(wavelengths):
        fwhm_scaled_px = psf_fwhm_func(wl) / arcsec_per_px_scaled
        convolved_model_scaled = convolve_image_with_gaussian(
            forward_model_scaled, fwhm_scaled_px
        )
        model_unscaled = convolved_model_scaled[idxs]
        psf_cube[wl_idx] = model_unscaled
    return psf_cube


# Default functions
def default_forward_model_func(body: planetmapper.BodyXY) -> np.ndarray:
    """
    Get forward model image for an observation.

    Illuminated pixels have a value of 1, unilluminated pixels or background sky have a
    value of 0.

    Args:
        body: planetmapper.BodyXY object to forward model.

    Returns:
        Numpy array containing forward model image.
    """
    return np.asarray(
        (body.get_emission_angle_img() < 90) & (body.get_incidence_angle_img() < 90),
        dtype=float,
    )


def default_psf_fwhm_func(wavelength: float) -> float:
    """
    PSF FWHM function from Law+2023.

    Gaussian PSF with FWHM given by θ = 0.033 (λ/micron) + 0.106 arcsec (Law+23).

    Args:
        wavelength: Wavelength in microns.

    Returns:
        PSF FWHM in arcseconds.
    """
    return 0.033 * wavelength + 0.106


# Helper functions
def create_scaled_body(
    observation: planetmapper.Observation,
    rescale: int,
    border: int,
) -> tuple[planetmapper.BodyXY, CoordinateTransform]:
    """
    Create a scaled version of the input observation to oversample the pixel grid when
    forward modelling.

    Args:
        observation: Observation to scale.
        rescale: Scaling factor to resize the image by.
        border: Number of pixels to add to the border of the image after rescaling.

    Returns:
        `(body, transform)` tuple where `body` is the scaled version of the input
        observation and `transform` is a function that can be used to transform
        pixel coordinates from the input observation to the scaled version.
    """
    body = planetmapper.BodyXY.from_body(observation.to_body())

    nx, ny = observation.get_img_size()
    body.set_img_size(
        nx * rescale + border * 2,
        ny * rescale + border * 2,
    )
    body.set_r0(observation.get_r0() * rescale)
    body.set_x0(observation.get_x0() * rescale + border)
    body.set_y0(observation.get_y0() * rescale + border)
    body.set_rotation(observation.get_rotation())

    transform: CoordinateTransform = lambda x: x * rescale + border
    return body, transform


def create_dynamically_scaled_body(
    observation: planetmapper.Observation,
    rescale: int | None = None,
    min_rescale: int = 5,
    min_disc_diameter_scaled: int = 100,
    buffer_unscaled: int = 6,
    max_target_diameter_to_skip: int = 1,
) -> tuple[planetmapper.BodyXY, CoordinateTransform, dict[str, int]]:
    """
    Automatically determine the size of the scaled observation to use for forward
    modelling.

    Args:
        observation: Observation to scale.
        rescale: Scaling factor to resize the image by. If `None`, this will be
            determined automatically using the following parametets.
        min_rescale: Minimum scaling factor to resize the image by.
        min_disc_diameter_scaled: Minimum pixel diameter of the disc the scaled image.
        buffer_unscaled: Number of pixels to add to the border of the image before
            rescaling.
        max_target_diameter_to_skip: If the target diameter is below this threshold,
            raise a `TargetTooSmallError`. This is to avoid scaling observations that
            are too small to have a useful PSF correction (these observations would also
            have a scaled observation that is very very large).

    Returns:
        `(body, transform, parameters)` tuple where `body` is the scaled version of the
        input observation, `transform` is a function that can be used to transform
        pixel coordinates from the input observation to the scaled version, and
        `parameters` is a dictionary containing the parameters used to scale the
        observation.

    Raises:
        TargetTooSmallError: If the target is too small to be scaled (as determoined by
            `max_target_diameter_to_skip`).
    """
    target_diameter_px = observation.get_r0() * 2 / observation.get_plate_scale_arcsec()
    if target_diameter_px < max_target_diameter_to_skip:
        raise TargetTooSmallError(
            f'Target diameter {target_diameter_px:.2f} px < '
            f'{max_target_diameter_to_skip} px'
        )

    if rescale is None:
        rescale = max(
            min_rescale,
            int(np.ceil(min_disc_diameter_scaled / target_diameter_px)),
        )
    buffer = buffer_unscaled * rescale
    body, transform = create_scaled_body(observation, rescale, buffer)
    parameters = {
        'rescale': rescale,
        'buffer': buffer,
        'min_rescale': min_rescale,
        'min_disc_diameter_scaled': min_disc_diameter_scaled,
        'buffer_unscaled': buffer_unscaled,
        'max_target_diameter_to_skip': max_target_diameter_to_skip,
    }
    return body, transform, parameters


def convolve_image_with_gaussian(image: np.ndarray, fwhm_px: float) -> np.ndarray:
    """
    Convolve an image with a Gaussian PSF.

    Args:
        image: Image to convolve.
        fwhm_px: FWHM of the Gaussian PSF in pixels.

    Returns:
        Convolved image.
    """
    sigma_px = fwhm_px / 2.35482  # sigma = FWHM / (2 * sqrt(2 * ln(2)))
    return gaussian_filter(image, sigma=sigma_px, mode='nearest')


def check_scaling(
    unscaled: planetmapper.BodyXY,
    scaled: planetmapper.BodyXY,
    transform: CoordinateTransform,
) -> None:
    """
    Run check to ensure unscaled and scaled bodies are consistent.
    """
    xy_coords = [
        (0, 0),
        (1, 3),
        (4.56, 7.89),
        (-0.5, 4.2),
    ]
    for x, y in xy_coords:
        radec_unscaled = unscaled.xy2radec(x, y)
        radec_scaled = scaled.xy2radec(transform(x), transform(y))
        if not np.allclose(radec_unscaled, radec_scaled):
            raise AssertionError(
                f'Failed consistency check for ({x}, {y}): '
                f'unscaled: {radec_unscaled}, scaled: {radec_scaled}'
            )
