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
__version__ = '0.1.0'
import argparse
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
        iterator = tqdm.tqdm(input_paths, desc=f'Defringing files')
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
    check_if_same_file: bool = True,
    print_info: bool = False,
) -> None:
    """
    Apply JWST pipeline 1D defringing to a single file.

    Uses `defringe_cube` to correct the SCI extension of the input file and saves the
    result to the output file.

    Args:
        input_path: Path to the input FITS file.
        output_path: Path to the output FITS file.
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
        )
        hdul['SCI'].data = cube_corrected  # type: ignore

        corrected_hdr = fits.Header()
        corrected_hdr.add_comment('Spaxels corrected by 1D defringing')
        corrected_hdr.add_comment('True values were corrected, False values were not')
        hdul.append(
            fits.ImageHDU(
                data=corrected_spaxels,
                header=corrected_hdr,
                name='1DFRINGE_CORRECTED',
            )
        )

        tools.add_header_reduction_note(hdul, '1D defringed')
        header['HIERARCH 1DFRINGE VERSION'] = (__version__, 'Software version')
        header['HIERARCH 1DFRINGE CHANNEL'] = (channel, 'Channel number')

        tools.check_path(output_path)
        hdul.writeto(output_path, overwrite=True)

    if print_info:
        print('  Done')


def defringe_cube(
    wavelengths: np.ndarray,
    cube: np.ndarray,
    channel: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply JWST pipeline 1D defringing to each spaxel of a 3D cube independently.

    If a spaxel cannot be corrected (if an erorr is raised while calling the pipeline's
    `fit_residual_fringes_1d` function), it will be left unchanged. The returned
    `corrected_spaxels` array can be used to identify which spaxels were corrected.

    Args:
        wavelengths: Wavelengths.
        cube: Data cube.
        channel: MIRI channel number.

    Returns:
        `(output_cube, corrected_spaxels)` tuple, where `output_cube` is the defringed
        data cube and `corrected_spaxels` is a boolean array indicating which spaxels
        were corrected (True values were corrected, False values were not).
    """
    output = cube.copy()
    corrected_spaxels = np.zeros(cube.shape[1:], dtype=bool)
    for i1 in range(cube.shape[1]):
        for i2 in range(cube.shape[2]):
            spectrum = cube[:, i1, i2]
            mask = np.isfinite(spectrum)
            try:
                output[mask, i1, i2] = fit_residual_fringes_1d(
                    flux=spectrum[mask],
                    wavelengths=wavelengths[mask],
                    channel=channel,
                )
                corrected_spaxels[i1, i2] = True
            except np.linalg.LinAlgError:
                pass
    return output, corrected_spaxels


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        argument_default=argparse.SUPPRESS,
    )
    parser.add_argument(
        'input-paths',
        nargs='+',
        help='Paths to the input FITS files.',
    )
    parser.add_argument(
        '--directory-out',
        '-d',
        help='Path to the output directory. Cannot be used in combination with output-path.',
    )
    parser.add_argument(
        '--output-path',
        '-o',
        help='Path to the output FITS file. Can only be used with a single input file.',
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
    if 'output_path' in args:
        if 'directory_out' in args:
            parser.error('Cannot use both output-path and directory-out')
        if len(args.input_paths) > 1:
            parser.error('Cannot use output-path with multiple input files')
        defringe_file(
            input_path=args.input_paths[0],
            output_path=args.output_path,
            check_if_same_file=args.check_if_same_file,
            print_info=args.print_info,
        )
    else:
        defringe_multiple(
            *args.input_paths,
            output_directory=args.directory_out,
            check_if_same_file=args.check_if_same_file,
            print_info=args.print_info,
        )


if __name__ == '__main__':
    main()
