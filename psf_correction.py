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


# Note: for MIRI, PSF FWHM is ~2.5 - ~3.5px for the longest wavelength of each channel


def correct_file(
    path: str,
    forward_model_func: ForwardModelFunc | None = None,
    psf_fwhm_func: FwhmFunc | None = None,
    mask_off_disc: bool = True,
    **kwargs,
) -> tuple[np.ndarray, np.ndarray]:
    """TODO"""
    if forward_model_func is None:
        forward_model_func = default_forward_model_func
    if psf_fwhm_func is None:
        psf_fwhm_func = default_psf_fwhm_func

    observation = planetmapper.Observation(path)
    try:
        body, transform, parameters = create_dynamically_scaled_body(
            observation, **kwargs
        )
    except TargetTooSmallError:
        return  # TODO

    with fits.open(path) as hdul:
        data_header: fits.Header = hdul['SCI'].header  # type: ignore
        cube: np.ndarray = hdul['SCI'].data  # type: ignore
        wavelengths = tools.get_wavelengths(data_header)

        psf_model = calculate_psf_modelled_cube(
            cube,
            wavelengths,
            body,
            transform,
            parameters['rescale'],
            forward_model_func,
            psf_fwhm_func,
        )

        cube_corrected = cube / psf_model
        forward_model = forward_model_func(observation)
        if mask_off_disc:
            cube_corrected[:, forward_model == 0] = np.nan

    return psf_model, cube_corrected


def calculate_psf_modelled_cube(
    cube: np.ndarray,
    wavelengths: np.ndarray,
    body: planetmapper.BodyXY,
    transform: CoordinateTransform,
    rescale: int,
    forward_model_func: ForwardModelFunc,
    psf_fwhm_func: FwhmFunc,
) -> np.ndarray:
    """TODO"""
    forward_model_scaled = forward_model_func(body)
    arcsec_per_px = body.get_plate_scale_arcsec()

    idxs1 = [transform(i) for i in range(cube.shape[1])]
    idxs2 = [transform(i) for i in range(cube.shape[2])]
    idxs = np.ix_(idxs1, idxs2)

    psf_cube = np.full_like(cube, np.nan)
    for wl_idx, wl in enumerate(wavelengths):
        fwhm_unscaled_px = psf_fwhm_func(wl) / arcsec_per_px
        convolved_model_scaled = convolve_image_with_gaussian(
            forward_model_scaled, fwhm_unscaled_px * rescale
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
    resize: int,
    border: int,
) -> tuple[planetmapper.BodyXY, CoordinateTransform]:
    """
    Create a scaled version of the input observation to oversample the pixel grid when
    forward modelling.

    Args:
        observation: Observation to scale.
        resize: Scaling factor to resize the image by.
        border: Number of pixels to add to the border of the image after rescaling.

    Returns:
        `(body, transform)` tuple where `body` is the scaled version of the input
        observation and `transform` is a function that can be used to transform
        pixel coordinates from the input observation to the scaled version.
    """
    body = planetmapper.BodyXY.from_body(observation.to_body())

    nx, ny = observation.get_img_size()
    body.set_img_size(
        nx * resize + border * 2,
        ny * resize + border * 2,
    )
    body.set_r0(observation.get_r0() * resize)
    body.set_x0(observation.get_x0() * resize + border)
    body.set_y0(observation.get_y0() * resize + border)
    body.set_rotation(observation.get_rotation())

    transform: CoordinateTransform = lambda x: x * resize + border
    return body, transform


def create_dynamically_scaled_body(
    observation: planetmapper.Observation,
    min_rescale: int = 2,
    min_disc_diameter_scaled: int = 100,
    buffer_unscaled: int = 10,
    max_target_diameter_to_skip: int = 1,
) -> tuple[planetmapper.BodyXY, CoordinateTransform, dict[str, int]]:
    """
    Automatically determine the size of the scaled observation to use for forward
    modelling.

    Args:
        observation: Observation to scale.
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


class TargetTooSmallError(Exception):
    """Raised when the target is too small to be scaled."""


def convolve_image_with_gaussian(image: np.ndarray, fwhm_px: float) -> np.ndarray:
    """TODO"""
    sigma_px = fwhm_px / 2.35482  # sigma = FWHM / (2 * sqrt(2 * ln(2)))
    return gaussian_filter(image, sigma=sigma_px, mode='nearest')
