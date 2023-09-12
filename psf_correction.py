import numpy as np
import planetmapper
from astropy.io import fits


def forward_model(path: str) -> np.ndarray:
    """
    Get forward model image for an observation.

    Illuminated pixels have a value of 1, unilluminated pixels or background sky have a
    value of 0.

    Args:
        path: Path to observation FITS file.

    Returns:
        Numpy array containing forward model image.
    """
    # TODO need to expand the image to include full planet rather than just the FOV
    # TODO probably want to oversample the image to model the PSF better
    observation = planetmapper.Observation(path)
    return np.asarray(
        (observation.get_emission_angle_img() < 90)
        & (observation.get_incidence_angle_img() < 90),
        dtype=float,
    )
