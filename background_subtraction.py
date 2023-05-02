#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from astropy.io import fits
import numpy as np
import tools
import os


def subtract_background(
    path_in: str,
    path_out: str,
    path_background: str,
) -> None:
    with fits.open(path_background) as hdul:
        background_cube: np.ndarray = hdul['SCI'].data  # type: ignore
        background_header: fits.Header = hdul['PRIMARY'].header  # type: ignore

    with fits.open(path_in) as hdul:
        cube = hdul['SCI'].data  # type: ignore
        header = hdul['PRIMARY'].header  #  type: ignore

        cube = cube - background_cube

        hdul['SCI'].data = cube  # type: ignore
        hdul['DQ'].data[np.isnan(hdul['SCI'].data)] += 1  #  type: ignore

        tools.add_header_reduction_note(hdul, 'Background subtracted')
        header['HIERARCH BACKGROUND'] = (True, 'Background subtracted')
        header['HIERARCH BACKGROUND OBS_ID'] = background_header['OBS_ID']
        bg_dir = os.path.split(path_background)[0]
        if len(bg_dir) > 48:
            bg_dir = '...' + bg_dir[-48:]
        header['HIERARCH BACKGROUND DIR'] = bg_dir
        bg_file = os.path.split(path_background)[1]
        if len(bg_file) > 47:
            bg_file = '...' + bg_file[-47:]
        header['HIERARCH BACKGROUND FILE'] = bg_file

        tools.check_path(path_out)
        hdul.writeto(path_out, overwrite=True)
