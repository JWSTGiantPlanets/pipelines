import sys
import tools
import numpy as np
from astropy.io import fits
import save_spx
import tqdm
import scipy.interpolate
from typing import Callable

DATA_KEYS = ['SCI', 'ERR']

BACKPLANE_KEYS = [
    'LAT',
    'LON',
    'LAT_LIMB',
    'LON_LIMB',
    'LAT_GRAPHIC',
    'PHASE',
    'EMISSION',
    'INCIDENCE',
    'AZIMUTH',
    'LOCALTIME',
    'DISTANCE_LIMB',
    'DISTANCE_RINGS',
    'LON_RINGS',
    'RA',
    'DEC',
    'RADIAL_VELOCITY',
    'DOPPLER',
]

HEADER_KEYS = ['MJD-AVG', 'PIXAR_A2', 'PIXAR_SR']


def generate_spx(path_in: str, path_out_template: str):
    with fits.open(path_in) as hdul:
        header = hdul['SCI'].header
        wavelengths = tools.get_wavelengths(header, axis=2)
        values = tools.get_wavelengths(header, axis=1, warn=False)
        for idx, val in enumerate(values):
            spectrum = hdul['SCI'].data[:, idx]
            error = hdul['ERR'].data[:, idx]

            if all(np.isnan(spectrum)):
                continue

            wl = wavelengths[np.isfinite(spectrum)]
            error = error[np.isfinite(spectrum)]
            spectrum = spectrum[np.isfinite(spectrum)]

            path_out = path_out_template.format(val=format(val, '+g'))
            save_spx.save_spx(
                path_out,
                wavelengths=wl,
                spectrum=spectrum,
                error=error,
                lon=hdul['LON'].data[idx],
                lat=hdul['LAT_GRAPHIC'].data[idx],
                phase=hdul['PHASE'].data[idx],
                emission=hdul['EMISSION'].data[idx],
                azimuth=hdul['AZIMUTH'].data[idx],
                pixelarea_arcsecsq=hdul['PIXAR_A2'].data[idx],
                pixelarea_steradians=hdul['PIXAR_SR'].data[idx],
            )


def average_grouped_data(
    path_in: str,
    path_out: str,
    reject_fraction: float = 1 / 3,
):
    with fits.open(path_in) as hdul:
        spectral_data = hdul['SCI'].data  # type: ignore
        mjd = hdul['MJD-AVG'].data  #  type: ignore
        masks = []
        for idx in range(spectral_data.shape[1]):
            mjd_mask = np.isfinite(mjd[idx, :])
            spectra = spectral_data[:, idx, mjd_mask]
            if spectra.shape[1] == 0:
                masks.append(None)
                continue
            avg = np.nanmedian(spectra, axis=1)
            avg = np.nan_to_num(avg)
            rms = np.array(
                [
                    tools.rms(avg, np.nan_to_num(spectra[:, i]))
                    for i in range(spectra.shape[1])
                ]
            )
            threshold = np.quantile(rms, 1 - reject_fraction)
            rms_mask = rms < threshold
            masks.append((mjd_mask, rms_mask))

        for k in DATA_KEYS + BACKPLANE_KEYS + HEADER_KEYS:
            data = hdul[k].data  #  type: ignore
            new = np.full(data.shape[:-1], np.nan)
            for idx, m in enumerate(masks):
                if m is None:
                    continue
                mjd_mask, rms_mask = m
                value_data = data[..., idx, :][..., mjd_mask][..., rms_mask]
                avg_data = np.nanmean(value_data, axis=-1)
                new[..., idx] = avg_data
            hdul[k].data = new  # type: ignore

            header = hdul[k].header  # type: ignore
            for n in [2, 3]:
                for prefix in ['CRPIX', 'CRVAL', 'CDELT', 'CTYPE', 'CUNIT']:
                    try:
                        header[f'{prefix}{n-1}'] = header[f'{prefix}{n}']
                        del header[f'{prefix}{n}']
                    except KeyError:
                        pass

        new = np.zeros(mjd.shape[:-1])
        for idx, m in enumerate(masks):
            if m is None:
                continue
            mjd_mask, rms_mask = m
            new[idx] = sum(rms_mask)

        hdul.append(
            fits.ImageHDU(
                data=new,
                header=hdul['MJD-AVG'].header,
                name='NPTS',
            )
        )
        hdul.writeto(path_out, overwrite=True)


def save_grouped_data(
    values: np.ndarray | tuple[float, ...],
    data: list[dict],
    wavelengths: np.ndarray,
    path: str,
    ctype: str,
):
    max_points = max(len(d['SCI']) for d in data)

    hdul = fits.HDUList([fits.PrimaryHDU()])
    for k in DATA_KEYS:
        cube = np.full((len(wavelengths), len(values), max_points), np.nan)
        for val_idx, d in enumerate(data):
            value_data = d[k]
            if len(value_data):
                cube[:, val_idx, : len(value_data)] = np.array(value_data).T

        header = fits.Header()
        header['CRPIX2'] = 1
        header['CRVAL2'] = values[0]
        header['CDELT2'] = values[1] - values[0]
        header['CTYPE2'] = ctype
        header['CUNIT2'] = 'deg'

        header['CRPIX3'] = 1
        header['CRVAL3'] = wavelengths[0]
        header['CDELT3'] = wavelengths[1] - wavelengths[0]
        header['CTYPE3'] = 'WAVE'
        header['CUNIT3'] = 'um'
        hdul.append(fits.ImageHDU(data=cube, header=header, name=k))

    for k in BACKPLANE_KEYS + HEADER_KEYS:
        cube = np.full((len(values), max_points), np.nan)
        for val_idx, d in enumerate(data):
            value_data = d[k]
            if len(value_data):
                cube[val_idx, : len(value_data)] = value_data
        header = fits.Header()
        header['CRPIX2'] = 1
        header['CRVAL2'] = values[0]
        header['CDELT2'] = values[1] - values[0]
        header['CTYPE2'] = ctype
        header['CUNIT2'] = 'deg'

        hdul.append(fits.ImageHDU(data=cube, header=header, name=k))
    hdul.writeto(path, overwrite=True)


def group_data(
    paths: list[str],
    values: np.ndarray | tuple[float, ...],
    value_key: str,
    dq_threshold=0.25,
    mask_condition: Callable[[fits.HDUList], np.ndarray] | None = None,
) -> tuple[list[dict], np.ndarray]:
    channel_data = None
    wavelengths = None
    data = [{k: [] for k in DATA_KEYS + BACKPLANE_KEYS + HEADER_KEYS} for v in values]
    for p in tqdm.tqdm(paths, desc='Reading data', leave=False):
        with fits.open(p) as hdul:
            header = hdul['PRIMARY'].header  # type: ignore
            data_header = hdul['SCI'].header  # type: ignore

            if channel_data is None:
                channel_data = miri_channel_data_from_header(header)
            if wavelengths is None:
                wavelengths = tools.get_wavelengths(data_header)

            assert channel_data == miri_channel_data_from_header(header)
            assert np.array_equal(wavelengths, tools.get_wavelengths(data_header))

            dq = hdul['DQ'].data  # type: ignore
            mask = np.where(dq != 0)

            data_cubes = {}
            for k in DATA_KEYS:
                data_cubes[k] = hdul[k].data  # type: ignore
                data_cubes[k][mask] = np.nan

            if mask_condition is not None:
                good_mask = mask_condition(hdul)
            else:
                good_mask = None

            for i1 in range(dq.shape[1]):
                for i2 in range(dq.shape[2]):
                    if good_mask is not None:
                        if not good_mask[i1, i2]:
                            continue
                    bad_frac = sum(dq[:, i1, i2] != 0) / dq.shape[0]
                    if bad_frac > dq_threshold:
                        continue

                    v = hdul[value_key].data[i1, i2]  # type: ignore
                    if np.isnan(v):
                        continue

                    idx = tools.nearest_idx(values, v)

                    for k in DATA_KEYS:
                        data[idx][k].append(data_cubes[k][:, i1, i2])

                    for k in BACKPLANE_KEYS:
                        img = hdul[k].data  # type: ignore
                        data[idx][k].append(img[i1, i2])

                    for k in HEADER_KEYS:
                        data[idx][k].append(data_header[k])

    assert wavelengths is not None
    return data, wavelengths


def subtract_zonal_average(
    cube: np.ndarray,
    lat_img: np.ndarray,
    zonal_averages: np.ndarray,
    hdr_zonal: fits.Header,
    interpolate: bool = True,
    fill_value=np.nan,
) -> np.ndarray:
    assert cube.shape[0] == zonal_averages.shape[0], 'Inconsistent wavelength lengths'
    latitudes = tools.get_wavelengths(hdr_zonal, warn=False, axis=1)

    interpolator = scipy.interpolate.interp1d(
        latitudes, zonal_averages, bounds_error=False, fill_value=fill_value
    )
    cube = cube.copy()
    for i1 in range(cube.shape[1]):
        for i2 in range(cube.shape[2]):
            lat = lat_img[i1, i2]
            if interpolate:
                cube[:, i1, i2] -= interpolator(lat)
            else:
                cube[:, i1, i2] -= zonal_averages[:, tools.nearest_idx(latitudes, lat)]
    return cube


def miri_channel_data_from_header(
    header: fits.Header,
) -> tuple[str, str]:
    """
    Get the channel and band from a FITS header.
    """
    return header['CHANNEL'].casefold(), header['BAND'].casefold()
