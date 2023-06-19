import numpy as np
from scipy.ndimage import generic_filter


def exp_despike(
    img, size: int = 3, nsigma: float = 5, iterations: int = 1
) -> np.ndarray:
    if iterations <= 0:
        return img
    med_image = generic_filter(img, np.nanmedian, size, mode='constant', cval=np.nan)
    diff_image = img - med_image
    bad_values = np.abs(diff_image) > np.nanstd(diff_image) * nsigma
    img[bad_values] = np.nan
    if iterations > 1 and np.any(bad_values):
        return exp_despike(img, nsigma=nsigma, iterations=iterations - 1)
    return img


def exp_despike_global(img, nsigma: float = 5, iterations: int = 1) -> np.ndarray:
    if iterations <= 0:
        return img
    bad_values = np.abs(img - np.nanmedian(img)) > np.nanstd(img) * nsigma
    img[bad_values] = np.nan
    if iterations > 1 and np.any(bad_values):
        return exp_despike_global(img, nsigma=nsigma, iterations=iterations - 1)
    return img


def sigma_filter(
    img: np.ndarray,
    size: int = 5,
    nsigma: float = 3,
    threshold: float = 0,
    fix_values: bool = False,
    iterations: int = 1,
) -> np.ndarray:
    if iterations <= 0:
        return img
    med_img = generic_filter(img, np.nanmedian, size, mode='nearest')
    diff_img = img - med_img
    std_img = generic_filter(diff_img, np.nanstd, size, mode='nearest')
    sigma_img = np.abs(diff_img) / std_img
    bad_values = np.logical_and(sigma_img > nsigma, np.abs(diff_img) > threshold)
    if fix_values:
        img[bad_values] = med_img[bad_values]
    else:
        img[bad_values] = np.nan
    if iterations > 1 and np.any(bad_values):
        # Only need to keep iterating if there were any bad values
        return sigma_filter(
            img,
            size=size,
            nsigma=nsigma,
            threshold=threshold,
            fix_values=fix_values,
            iterations=iterations - 1,
        )
    return img
