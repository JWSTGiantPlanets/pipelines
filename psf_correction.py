from typing import Callable, TypeAlias

import numpy as np
import planetmapper
from astropy.io import fits

# Transform from original to scaled pixel coordinates:
CoordinateTransform: TypeAlias = Callable[[int], int]


def forward_model(body: planetmapper.BodyXY) -> np.ndarray:
    """
    Get forward model image for an observation.

    Illuminated pixels have a value of 1, unilluminated pixels or background sky have a
    value of 0.

    Args:
        body: planetmapper.BodyXY object to forward model.

    Returns:
        Numpy array containing forward model image.
    """
    return np.asarray(
        (body.get_emission_angle_img() < 90) & (body.get_incidence_angle_img() < 90),
        dtype=float,
    )


def create_scaled_body(
    observation: planetmapper.Observation,
    resize: int,
    border: int,
) -> tuple[planetmapper.BodyXY, CoordinateTransform]:
    """
    Create a scaled version of the input observation to oversample the pixel grid when
    forward modelling.

    Args:
        observation: Observation to scale.
        resize: Scaling factor to resize the image by.
        border: Number of pixels to add to the border of the image after rescaling.

    Returns:
        `(body, transform)` tuple where `body` is the scaled version of the input
        observation and `transform` is a function that can be used to transform
        pixel coordinates from the input observation to the scaled version.
    """
    body = planetmapper.BodyXY.from_body(observation.to_body())

    nx, ny = observation.get_img_size()
    body.set_img_size(nx * resize + border * 2, ny * resize + border * 2)
    body.set_r0(observation.get_r0() * resize)
    body.set_x0(observation.get_x0() * resize + border)
    body.set_y0(observation.get_y0() * resize + border)
    body.set_rotation(observation.get_rotation())

    transform: CoordinateTransform = lambda x: x * resize + border
    return body, transform


def create_dynamically_scaled_body(
    observation: planetmapper.Observation,
) -> tuple[planetmapper.BodyXY, CoordinateTransform]:
    # TODO generate scaled body using appropriate values for resize & border
    pass
