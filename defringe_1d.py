#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Apply JWST pipeline 1D defringing to FITS files.

example usage:
    python3 defringe_1d.py data.fits -o data_defringed.fits
    python3 defringe_1d.py data1.fits data2.fits -d output_dir
    python3 defringe_1d.py data/input/*.fits -d data/output
    python3 defringe_1d.py data.fits -o data.fits --no-check-if-same-file
"""
__version__ = '1.1.1'
import argparse
from collections.abc import Collection
from pathlib import Path

import numpy as np
import tqdm
from astropy.io import fits
from jwst.residual_fringe.utils import fit_residual_fringes_1d

import tools


def defringe_multiple(
    *input_paths: str | Path,
    output_directory: str | Path,
    print_info: bool = False,
    **kwargs,
) -> None:
    """
    Apply JWST pipeline 1D defringing to multiple FITS files.

    Output files will be saved in `output_directory` and will have the same names as the
    input files.

    Args:
        input_paths: Paths to the input FITS files.
        output_directory: Path to the output directory.
        print_info: If True, print logging information.
        **kwargs: Additional keyword arguments to pass to `defringe_file`.
    """
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    if print_info:
        iterator = tqdm.tqdm(input_paths, desc='Defringing files')
    else:
        iterator = input_paths

    for path_in in iterator:
        path_in = Path(path_in)
        path_out = output_directory / path_in.name
        defringe_file(path_in, path_out, **kwargs)


def defringe_file(
    input_path: str | Path,
    output_path: str | Path,
    *,
    wavelength_ranges: Collection[tuple[float, float]] | None = None,
    check_if_same_file: bool = True,
    print_info: bool = False,
) -> None:
    """
    Apply JWST pipeline 1D defringing to a single file.

    Uses `defringe_cube` to correct the SCI extension of the input file and saves the
    result to the output file. The corrected spaxels are saved in a new extension called
    'FRINGE1D_CORRECTED'.

    Args:
        input_path: Path to the input FITS file.
        output_path: Path to the output FITS file.
        wavelength_ranges: List of wavelength ranges to use for defringing. See the
            `defringe_spectrum` function for more details.
        check_if_same_file: If True, raise a ValueError if `input_path` and
            `output_path` are the same.
        print_info: If True, print logging information.
    """
    input_path = Path(input_path).resolve()
    output_path = Path(output_path).resolve()
    if check_if_same_file and input_path == output_path:
        raise ValueError('Input and output paths are the same.')
    if print_info:
        print(f'Correcting {input_path} -> {output_path}')

    with fits.open(input_path) as hdul:
        header = hdul['PRIMARY'].header  #  type: ignore
        cube = hdul['SCI'].data  # type: ignore
        data_header = hdul['SCI'].header  #  type: ignore

        channel = int(header['CHANNEL'])
        wavelengths = tools.get_wavelengths(data_header)

        cube_corrected, corrected_spaxels = defringe_cube(
            wavelengths=wavelengths,
            cube=cube,
            channel=channel,
            wavelength_ranges=wavelength_ranges,
            print_info=print_info,
        )
        hdul['SCI'].data = cube_corrected  # type: ignore

        corrected_hdr = fits.Header()
        corrected_hdr.add_comment('Spaxels corrected by 1D defringing')
        corrected_hdr.add_comment('True values were corrected, False values were not')
        hdul.append(
            fits.ImageHDU(
                data=corrected_spaxels.astype(np.uint8),
                header=corrected_hdr,
                name='FRINGE1D_CORRECTED',
            )
        )

        tools.add_header_reduction_note(hdul, 'Residual fringe corrected (1D)')
        header['HIERARCH FRINGE1D VERSION'] = (__version__, 'Software version')
        header['HIERARCH FRINGE1D CHANNEL'] = (channel, 'Channel number')
        header['HIERARCH FRINGE1D USED_WAVL_RANGES'] = (
            wavelength_ranges is not None,
            'Limited defringing to specific wavelength ranges',
        )
        if wavelength_ranges is not None:
            header['HIERARCH FRINGE1D NUM_WAVL_RANGES'] = (
                len(wavelength_ranges),
                'Number of wavelength ranges used for defringing',
            )
            for i, (start, stop) in enumerate(wavelength_ranges):
                header[f'HIERARCH FRINGE1D WAVL_RANGE_{i}_START'] = start
                header[f'HIERARCH FRINGE1D WAVL_RANGE_{i}_STOP'] = stop

        tools.check_path(output_path)
        hdul.writeto(output_path, overwrite=True)


def defringe_cube(
    wavelengths: np.ndarray,
    cube: np.ndarray,
    channel: int,
    *,
    wavelength_ranges: Collection[tuple[float, float]] | None = None,
    print_info: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply JWST pipeline 1D defringing to each spaxel of a 3D cube independently.

    If a spaxel cannot be corrected (if an error is raised while calling the pipeline's
    `fit_residual_fringes_1d` function), it will be left unchanged. The returned
    `corrected_spaxels` array can be used to identify which spaxels were corrected.

    Args:
        wavelengths: Wavelengths.
        cube: Data cube.
        channel: MIRI channel number.
        wavelength_ranges: List of wavelength ranges to use for defringing. See the
            `defringe_spectrum` function for more details.
        print_info: If True, print logging information.

    Returns:
        `(output_cube, corrected_spaxels)` tuple, where `output_cube` is the defringed
        data cube and `corrected_spaxels` is a boolean array indicating which spaxels
        were corrected (True values were corrected, False values were not).
    """
    output = cube.copy()
    corrected_spaxels = np.zeros(cube.shape[1:], dtype=bool)

    indices = np.ndindex(cube.shape[1:])
    if print_info:
        indices = tqdm.tqdm(list(indices), desc='Defringing')

    for i1, i2 in indices:
        output[:, i1, i2], corrected_spaxels[i1, i2] = defringe_spectrum(
            wavelengths=wavelengths,
            spectrum=cube[:, i1, i2],
            channel=channel,
            wavelength_ranges=wavelength_ranges,
        )
    return output, corrected_spaxels


def defringe_spectrum(
    *,
    wavelengths: np.ndarray,
    spectrum: np.ndarray,
    channel: int,
    wavelength_ranges: Collection[tuple[float, float]] | None = None,
) -> tuple[np.ndarray, bool]:
    """
    Apply JWST pipeline 1D defringing to a 1D spectrum. This is a wrapper around the
    pipeline's `fit_residual_fringes_1d` function, with some extra checks and handling
    of errors.

    Args:
        wavelengths: Wavelengths.
        spectrum: 1D spectrum.
        channel: MIRI channel number.
        wavelength_ranges: List of wavelength ranges to use for defringing. This should
            be a list of tuples, where each tuple contains the start and end wavelengths
            of a range to use. For example, if `wavelength_ranges=[(5, 6), (7.5, 8)]`,
            the defringing will be performed using the data between 5 and 6 microns and
            between 7.5 and 8 microns (inclusive). Note that all provided wavelength
            ranges are defringed simultaneously. Any data outside of these ranges is not
            used as an input to the defringing algorithm, and will be left unchanged by
            this function. If None, the entire wavelength range will be used.

    Returns:
        `(output_spectrum, was_defringed)` tuple, where `output_spectrum` is the
        defringed spectrum and `was_defringed` is a boolean indicating whether the
        spectrum was defringed (True) or not (False).
    """
    output = spectrum.copy()
    was_defringed = False

    mask = np.isfinite(spectrum)
    if wavelength_ranges is not None:
        wavelength_range_mask = np.zeros_like(wavelengths, dtype=bool)
        for start, stop in wavelength_ranges:
            wavelength_range_mask |= (wavelengths >= start) & (wavelengths <= stop)
        mask &= wavelength_range_mask

    if any(mask):
        try:
            output[mask] = fit_residual_fringes_1d(
                flux=spectrum[mask],
                wavelength=wavelengths[mask],
                channel=channel,
            )
            was_defringed = True
        except (np.linalg.LinAlgError, ValueError, IndexError):
            pass
    return output, was_defringed


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        argument_default=argparse.SUPPRESS,
    )
    parser.add_argument(
        'input_paths',
        nargs='+',
        help='Paths to the input FITS files.',
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--directory-out',
        '-d',
        help='Path to the output directory. Cannot be used in combination with output-path.',
    )
    group.add_argument(
        '--output-path',
        '-o',
        help='Path to the output FITS file. Can only be used with a single input file.',
    )
    parser.add_argument(
        '--wavelength-range',
        '-w',
        nargs=2,
        type=float,
        action='append',
        default=None,
        metavar=('START', 'STOP'),
        help=(
            'Wavelength ranges to defringe, can be used multiple times. '
            'For example, use `-w 1 2 -w 3 4` to only defringe between 1 and 2 microns '
            'and between 3 and 4 microns. If not provided, the entire wavelength range '
            'will be used.'
        ),
    )
    parser.add_argument(
        '--print-info',
        action=argparse.BooleanOptionalAction,
        help='Print logging information.',
        default=True,
    )
    parser.add_argument(
        '--no-check-if-same-file',
        action='store_false',
        default=True,
        dest='check_if_same_file',
        help='Allow the input and output paths to be the same.',
    )
    args = parser.parse_args()
    kwargs = {k: v for k, v in vars(args).items() if k not in {'input_paths'}}
    if 'output_path' in args:
        if len(args.input_paths) > 1:
            parser.error('Cannot use output-path with multiple input files')
        defringe_file(input_path=args.input_paths[0], **kwargs)
    else:
        defringe_multiple(*args.input_paths, **kwargs)


if __name__ == '__main__':
    main()
