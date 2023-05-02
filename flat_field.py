#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from astropy.io import fits
import numpy as np
import tools

GENERATION_HEADER_PREFIX = 'FLAT GENERATION'

def apply_flat(
    path_in: str,
    path_out: str,
    path_flat: str,
) -> None:
    flat_hdus: list[tuple[np.ndarray, fits.Header]] = []
    with fits.open(path_flat) as hdul:
        flat_header: fits.Header = hdul[0].header  # type: ignore
        flat_cube_raw: np.ndarray = hdul[0].data  # type: ignore
        for hdu in hdul[1:]:
            flat_hdus.append((hdu.data, hdu.header))  #  type: ignore

    with fits.open(path_in) as hdul:
        cube = hdul['SCI'].data  # type: ignore
        header = hdul['PRIMARY'].header  #  type: ignore

        flat_cube = get_consistent_shaped_flat(cube, flat_cube_raw)
        cube = cube / flat_cube

        hdul['SCI'].data = cube  # type: ignore
        hdul['DQ'].data[np.isnan(hdul['SCI'].data)] += 1  #  type: ignore

        tools.add_header_reduction_note(hdul, 'Flat field corrected')
        for k, v in flat_header.items():
            if GENERATION_HEADER_PREFIX in k:
                if not k.upper().startswith('HIERARCH'):
                    k = f'HIERARCH {k}'
                header[k] = v
        header['HIERARCH FLAT RESHAPED'] = (
            flat_cube.shape != flat_cube_raw.shape,
            'Flat cube reshaped to match science cube',
        )
        header['HIERARCH FLAT ORIGINAL_SHAPE'] = str(flat_cube_raw.shape)

        for hdu in flat_hdus:
            hdul.append(fits.ImageHDU(data=hdu[0], header=hdu[1]))
        tools.check_path(path_out)
        hdul.writeto(path_out, overwrite=True)


def get_consistent_shaped_flat(cube: np.ndarray, flat: np.ndarray) -> np.ndarray:
    if cube.shape == flat.shape:
        return flat
    if cube.shape[0] != flat.shape[0]:
        raise ValueError(
            'Flat cube and science cube have different number of wavelengths'
        )
    for dim in [1, 2]:
        diff = cube.shape[dim] - flat.shape[dim]
        pad, remainder = divmod(diff, 2)
        if remainder:
            raise ValueError(
                '\n'.join(
                    [
                        f'Shape difference {diff} in dimension {dim} is not divisible by 2',
                        f'Cube shape: {cube.shape}',
                        f'Flat shape: {flat.shape}',
                    ]
                )
            )
        if pad > 0:
            padding = np.zeros((len(cube.shape), 2), dtype=int)
            padding[dim, :] = pad
            flat = np.pad(flat, padding, mode='constant', constant_values=np.nan)
        elif pad < 0:
            pad = abs(pad)
            slices = [slice(None)] * len(cube.shape)
            slices[dim] = slice(pad, -pad)
            flat = flat[tuple(slices)]
    return flat
