#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NIRSPEC pipeline

TODO
- Reduction stages 1-3
- Navigate
- Desaturation (optional)
- Despike
- Plot
- Animate
- Add note about chmodding all the touched files
"""
import argparse
import datetime
import glob
import json
import os
from typing import Any, Literal, TypeAlias

import tqdm

import background_subtraction
import desaturate_data
import despike_data
import flat_field
import jwst_summary_animation
import jwst_summary_plots
import navigate_jwst_observations
import reduce_jwst_miri
import remove_groups
import tools
from parallel_tools import runmany
from tools import KeepMissingDict


STEP: TypeAlias = Literal[
    'remove_groups',
    'reduce',
    'navigate',
    'desaturate',
    'despike',
    'background',
    'plot',
    'animate',
]
STEPS: list[STEP] = [
    'remove_groups',
    'reduce',
    'navigate',
    'desaturate',
    'despike',
    'background',
    'plot',
    'animate',
]


def run_pipeline(
    root_path: str,
    *,
    desaturate: bool = True,
    groups_to_use: list[int] | None = None,
    background_path: str | None = None,
    ndither: int = 4,
    parallel: float | bool = False,
    basic_navigation: bool = False,
    skip_steps: list[STEP] | set[STEP] | None = None,
    start_step: STEP | None = None,
    end_step: STEP | None = None,
) -> None:
    pass


def main(*args):
    pass
    # TODO argparse stuff


if __name__ == '__main__':
    main()
