import datetime
import math
import statistics
from typing import TypeAlias

import numpy as np
import tqdm
from astropy.io import fits

import flat_field
import tools

VARIABLE_HEADER_KEYS = ['DATE', 'DATASET']
HEADER_PREFIX = flat_field.GENERATION_HEADER_PREFIX

CorrespondingPixels: TypeAlias = dict[
    tuple[int, tuple[int, int]], set[tuple[int, tuple[int, int]]]
]
PixelRatios: TypeAlias = dict[tuple[int, int], list[tuple[tuple[int, int], float]]]


def construct_flat_cube(
    paths: list[str],
    lat_bin_size_factor: float,
    bin_aspect: float,
    emission_cutoff: float,
    iterations: int,
    nsigma: float,
    sigma_iterations: int,
    snr_threshold: float,
    max_correction: float,
    **tqdm_kw,
) -> np.ndarray:
    cubes: list[np.ndarray] = []
    error_cubes: list[np.ndarray] = []
    masks: list[np.ndarray] = []
    coord_images: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    shape = None
    for p in paths:
        with fits.open(p) as hdul:
            cube = hdul['SCI'].data  # type: ignore
            error_cube = hdul['ERR'].data  # type: ignore
            if shape is None:
                shape = cube.shape
            assert cube.shape == shape, 'Cube shape mismatch'
            try:
                lat_img = hdul['LAT_GRAPHIC'].data  # type: ignore
                lon_img = hdul['LON'].data  # type: ignore
            except KeyError:
                lat_img = hdul['LAT_PGR'].data  # type: ignore
                lon_img = hdul['LON_WEST'].data  # type: ignore
            emission_img = hdul['EMISSION'].data  # type: ignore
            try:
                dq = hdul['DQ'].data  #  type: ignore
                cube[dq != 0] = np.nan
            except KeyError:
                pass
            cube[cube == 0] = np.nan
            error_cube[np.isnan(cube)] = np.nan
            mask = emission_img < emission_cutoff

            error_cubes.append(error_cube)
            cubes.append(cube)
            masks.append(mask)
            coord_images.append((lon_img, lat_img, emission_img))
    assert shape is not None, 'No paths provided'

    flat_cube = np.ones(shape)

    corresponding_pixels = calculate_corresponding_pixels(
        coord_images, lat_bin_size_factor, bin_aspect, emission_cutoff
    )

    with tools.ignore_warnings('Mean of empty slice', 'All-NaN slice encountered'):
        for idx in tqdm.tqdm(range(shape[0]), **tqdm_kw):
            images = [cube[idx] for cube in cubes]
            erorrs = [error_cube[idx] for error_cube in error_cubes]
            signal = np.nanmean(
                [np.nanmean(img[mask]) for img, mask in zip(images, masks)]
            )
            noise = np.nanmean(
                [np.nanmean(err[mask]) for err, mask in zip(erorrs, masks)]
            )
            if signal / noise < snr_threshold:
                continue

            for _ in range(sigma_iterations):
                sigma = np.nanstd(images)
                median = np.nanmedian(images)
                for img in images:
                    img[img > median + nsigma * sigma] = np.nan
                    img[img < median - nsigma * sigma] = np.nan

            flat_img = calculate_flat_image(
                images, corresponding_pixels, iterations, max_correction
            )
            flat_cube[idx] = flat_img

    return flat_cube


def calculate_corresponding_pixels(
    coord_images: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    lat_bin_size_factor: float,
    bin_aspect: float,
    emission_cutoff: float,
) -> CorrespondingPixels:
    lat_images = [lat_img for lon_img, lat_img, emission_img in coord_images]
    lat_diff = (
        np.nanmedian(np.diff(lat_images, axis=1)) ** 2
        + np.nanmedian(np.diff(lat_images, axis=2)) ** 2
    ) ** 0.5
    lonlat_comparison_radius = lat_diff * lat_bin_size_factor

    corresponding_pixels: CorrespondingPixels = {}
    for img_idx, (lon_img, lat_img, emission_img) in enumerate(coord_images):
        for idxs, emission in np.ndenumerate(emission_img):
            if not emission < emission_cutoff:
                # Automatically checks for NaNs too
                continue
            lon = lon_img[idxs]
            lat = lat_img[idxs]
            key = (img_idx, idxs)
            corresponding_pixels.setdefault(key, set())
            for img2_idx, (lon_img2, lat_img2, emission_img2) in enumerate(
                coord_images
            ):
                for idxs2, emission2 in np.ndenumerate(emission_img2):
                    key2 = (img2_idx, idxs2)
                    if key2 in corresponding_pixels:
                        if key in corresponding_pixels[key2]:
                            corresponding_pixels[key].add(key2)
                    else:
                        if not emission2 < emission_cutoff:
                            continue
                        dlat = lat - lat_img2[idxs2]
                        # Do partial checks first for speed
                        if dlat > lonlat_comparison_radius:
                            continue
                        dlon = (lon - lon_img2[idxs2]) / bin_aspect
                        if dlon > lonlat_comparison_radius:
                            continue
                        if dlat**2 + dlon**2 < lonlat_comparison_radius**2:
                            corresponding_pixels[key].add(key2)
    return corresponding_pixels


def calculate_flat_image(
    images: list[np.ndarray],
    corresponding_pixels: CorrespondingPixels,
    iterations: int,
    max_correction: float,
) -> np.ndarray:
    flat = np.full(images[0].shape, np.nan)
    pixel_ratios = calculate_pixel_ratios(images, corresponding_pixels)

    # Automatically choose a starting pixel near the centre of the image which has data
    start_idxs = tuple(x // 2 for x in flat.shape)
    index_options = iter(
        sorted(
            np.ndindex(flat.shape),
            key=lambda x: float(np.linalg.norm(np.array(x) - np.array(start_idxs))),
        )
    )
    try:
        while start_idxs not in pixel_ratios:
            start_idxs = next(index_options)
    except StopIteration:
        # No valid start pixels found, so just fill flat with NaN
        return flat
    flat[start_idxs] = 1
    for iteration in range(iterations):
        flat_new = np.full_like(flat, np.nan)
        for idxs, ratios in pixel_ratios.items():
            # Faster than np.nanmedian
            values = [
                v * ratio for idxs2, ratio in ratios if not math.isnan(v := flat[idxs2])
            ]
            if values:
                value = statistics.median(values)
                if iteration < 1 or 1 / max_correction < value < max_correction:
                    # Allow the first step to proceeed without a check to avoid bad
                    # start pixel preventing flat from propagating
                    flat_new[idxs] = value
        flat = flat_new / np.nanmean(flat)
    return flat


def calculate_pixel_ratios(
    images: list[np.ndarray], corresponding_pixels: CorrespondingPixels
) -> PixelRatios:
    pixel_ratios: PixelRatios = {}
    for (img_idx, idxs), corresponding in corresponding_pixels.items():
        value = images[img_idx][idxs]
        if math.isnan(value) or value <= 0:
            continue
        pixel_ratios.setdefault(idxs, [])
        for img2_idx, idxs2 in corresponding:
            value2 = images[img2_idx][idxs2]
            if math.isnan(value2) or value2 <= 0:
                continue
            pixel_ratios[idxs].append((idxs2, value / value2))
    pixel_ratios = {idxs: r for idxs, r in pixel_ratios.items() if r}
    return pixel_ratios


def do_tile(
    paths: list[str],
    path_out: str,
    channel: str,
    band: str,
    fringe: str,
    dataset: str,
    *,
    lat_bin_size_factor: float = 1,
    bin_aspect: float = 1,
    emission_cutoff: float = 75,
    iterations: int = 50,
    nsigma: float = 5,
    sigma_iterations: int = 5,
    snr_threshold: float = 10,
    max_correction: float = 1.5,
):
    flat_cube = construct_flat_cube(
        paths,
        lat_bin_size_factor=lat_bin_size_factor,
        bin_aspect=bin_aspect,
        emission_cutoff=emission_cutoff,
        iterations=iterations,
        nsigma=nsigma,
        sigma_iterations=sigma_iterations,
        snr_threshold=snr_threshold,
        max_correction=max_correction,
        desc=dataset,
        leave=False,
    )

    header = fits.Header()
    header[f'HIERARCH {HEADER_PREFIX} DATE'] = datetime.datetime.now().isoformat()
    header[f'HIERARCH {HEADER_PREFIX} LAT_BIN_SIZE_FACTOR'] = lat_bin_size_factor
    header[f'HIERARCH {HEADER_PREFIX} BIN_ASPECT'] = bin_aspect
    header[f'HIERARCH {HEADER_PREFIX} EMISSION_CUTOFF'] = emission_cutoff
    header[f'HIERARCH {HEADER_PREFIX} ITERATIONS'] = iterations
    header[f'HIERARCH {HEADER_PREFIX} NSIGMA'] = nsigma
    header[f'HIERARCH {HEADER_PREFIX} SIGMA_ITERATIONS'] = sigma_iterations
    header[f'HIERARCH {HEADER_PREFIX} SNR_THRESHOLD'] = snr_threshold
    header[f'HIERARCH {HEADER_PREFIX} MAX_CORRECTION'] = max_correction
    header[f'HIERARCH {HEADER_PREFIX} DATASET'] = dataset
    header[f'HIERARCH {HEADER_PREFIX} CHANNEL'] = channel
    header[f'HIERARCH {HEADER_PREFIX} BAND'] = band
    header[f'HIERARCH {HEADER_PREFIX} DEFRINGED'] = bool(fringe)

    hdul = make_hdul(flat_cube, header)

    tools.check_path(path_out)
    hdul.writeto(path_out, overwrite=True)


def make_hdul(flat_cube: np.ndarray, header: fits.Header):
    corrected = np.any(flat_cube != 1, axis=(1, 2))
    header_corrected = fits.Header()
    header_corrected['DESC'] = 'Wavelengths corrected (1) or skipped (0) by flat field'

    flat_std = np.nanstd(flat_cube, axis=(1, 2))
    flat_std[~corrected] = np.nan
    header_std = fits.Header()
    header_std['DESC'] = 'Standard deviation of flat field'

    hdul = fits.HDUList(
        [
            fits.PrimaryHDU(data=flat_cube, header=header),
            fits.ImageHDU(
                data=np.asarray(corrected, dtype=float),
                header=header_corrected,
                name='FLAT_CORRECTED',
            ),
            fits.ImageHDU(data=flat_std, header=header_std, name='FLAT_STD'),
        ]
    )

    return hdul


def merge_flats(paths: list[str], path_out: str, *, merge_nan_threshold: float = 0.5):
    data: list[tuple[np.ndarray, fits.Header]] = [
        fits.getdata(p, header=True) for p in paths
    ]  # type: ignore
    cubes = [d[0] for d in data]
    headers = [d[1] for d in data]

    variable_header_keys = [
        f'HIERARCH {HEADER_PREFIX} {k}' for k in VARIABLE_HEADER_KEYS
    ]

    standardise = lambda s: s.upper().replace('HIERARCH ', '')
    for h in headers:
        for k, v in h.items():
            k = k.upper()
            if HEADER_PREFIX in k and standardise(k) not in [
                standardise(s) for s in variable_header_keys
            ]:
                assert v == headers[0][k], f'{k} not consistent'

    with tools.ignore_warnings('Mean of empty slice'):
        flat_cube = np.full_like(cubes[0], np.nan)
        for idx, _ in enumerate(flat_cube):
            flat_cube[idx] = merge_images(
                [cube[idx] for cube in cubes], merge_nan_threshold=merge_nan_threshold
            )

    header = headers[0].copy()

    for k in variable_header_keys:
        for idx, h in enumerate(headers):
            header[f'{k}_{idx+1}'] = h[k]
        try:
            del header[k]
        except KeyError:
            pass

    header[f'HIERARCH {HEADER_PREFIX} MERGE_NAN_THRESHOLD'] = merge_nan_threshold
    header[f'HIERARCH {HEADER_PREFIX} N_MERGED'] = len(cubes)
    header[f'HIERARCH {HEADER_PREFIX} MERGE DATE'] = datetime.datetime.now().isoformat()

    hdul = make_hdul(flat_cube, header)

    tools.check_path(path_out)
    hdul.writeto(path_out, overwrite=True)


def merge_images(images: list[np.ndarray], merge_nan_threshold: float) -> np.ndarray:
    mask = np.all(np.isfinite(images), axis=0)
    masked_images = [img[mask] for img in images]
    avg = np.nanmedian(masked_images)
    ratios = [np.nanmedian(img) / avg for img in masked_images]
    scaled = [img / ratios[idx] for idx, img in enumerate(images)]
    merged = np.nanmedian(scaled, axis=0)
    merged = merged / np.nanmean(merged)

    nan_frac = np.sum(np.isnan(merged)) / merged.size
    if nan_frac > merge_nan_threshold:
        merged[:] = 1
    return merged
