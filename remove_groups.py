import os

import numpy as np
from astropy.io import fits

import tools


def remove_groups_from_file(path: str, groups_to_use: list[int] | None = None) -> None:
    with fits.open(path) as hdul:
        ngroups_full = hdul[0].header['NGROUPS']  #  type: ignore
        data = hdul['SCI'].data  #  type: ignore
        original_header = hdul[0].header.copy()  #  type: ignore
        for ngroups in range(1, ngroups_full):
            if groups_to_use is not None and ngroups not in groups_to_use:
                continue
            shape = (data.shape[0], ngroups, data.shape[2], data.shape[3])
            reduced_data = np.zeros(shape, dtype=data.dtype)
            for idx in range(ngroups):
                reduced_data[:, idx, :, :] = data[:, idx, :, :]
            reduced_data = reduced_data.astype('uint16')

            hdul['SCI'].data = reduced_data  #  type: ignore

            hdul[0].header = original_header.copy()  #  type: ignore
            hdul[0].header['NGROUPS'] = ngroups  #  type: ignore
            tools.add_header_reduction_note(
                hdul, f'Using {ngroups}/{ngroups_full} groups'
            )
            root, filename = os.path.split(path)
            root, stage0 = os.path.split(root)
            path_out = os.path.join(
                root, 'groups', f'{ngroups}_groups', stage0, filename
            )
            tools.check_path(path_out)
            hdul.writeto(path_out, overwrite=True)
