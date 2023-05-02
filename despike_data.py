import tools
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import tqdm
import image_filters


def despike_cube(
    path: str,
    path_out: str,
    size: int = 9,
    nsigma: float = 7,
    iterations: int = 3,
    progress_bar: bool = True,
):
    with fits.open(path) as hdul:
        cube: np.ndarray = hdul['SCI'].data  # type: ignore
        dq = hdul['DQ'].data  # type: ignore
        cube[dq != 0] = np.nan
        with tools.ignore_warnings('All-NaN slice encountered'):
            cube_to_loop = cube
            if progress_bar:
                cube_to_loop = tqdm.tqdm(cube_to_loop, 'Despiking cube', leave=False)

            for idx, img in enumerate(cube_to_loop):
                cube[idx] = image_filters.exp_despike(
                    img,
                    size=size,
                    nsigma=nsigma,
                    iterations=iterations,
                )

        hdul['SCI'].data = cube  # type: ignore
        hdul['DQ'].data[np.isnan(hdul['SCI'].data)] += 1  #  type: ignore
        tools.add_header_reduction_note(hdul, 'Despiked')
        header = hdul['PRIMARY'].header  # type: ignore

        header['HIERARCH CLEAN'] = (True, 'Set bad pixels to NaN')
        header['HIERARCH CLEAN METHOD'] = 'exp_despike'
        header['HIERARCH CLEAN SIZE'] = (size, 'Averaging window size')
        header['HIERARCH CLEAN NSIGMA'] = (nsigma, 'Standard deviation threshold')
        header['HIERARCH CLEAN ITERATIONS'] = iterations

        tools.check_path(path_out)
        hdul.writeto(path_out, overwrite=True)
