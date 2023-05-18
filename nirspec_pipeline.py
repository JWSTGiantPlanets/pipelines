#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NIRSPEC pipeline

TODO
- Remove groups (optional)
- Reduction stages 3
- Navigate
- Desaturation (optional)
- Despike
- Plot
- Animate
- Add note about chmodding all the touched files
"""
import glob
import os
from typing import Any, Literal, TypeAlias

from jwst.pipeline import Detector1Pipeline, Spec2Pipeline

from parallel_tools import runmany
from tools import check_path, log

# Central record of the filepaths for the various steps of the reduction pipeline


STEP: TypeAlias = Literal[
    'remove_groups',
    'stage1',
    'stage2',
    'stage3',
    'navigate',
    'desaturate',
    'despike',
    'plot',
    'animate',
]
STEPS: list[STEP] = [
    'remove_groups',
    'stage1',
    'stage2',
    'stage3',
    'navigate',
    'desaturate',
    'despike',
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
    stage1_kwargs: dict[str, Any] | None = None,
    stage2_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    # Process args
    parallel_kwargs = dict(
        parallel_frac=parallel,
        timeout=60 * 60,
    ) | (parallel_kwargs or {})
    reduction_parallel_kwargs = (
        parallel_kwargs
        | dict(
            start_delay=10,
            parallel_job_kw=dict(
                caught_error_wait_time=15,
                caught_error_wait_time_frac=1,
                caught_error_wait_time_max=600,
            ),
        )
        | (reduction_parallel_kwargs or {})
    )

    if skip_steps is None:
        skip_steps = []
    for s in skip_steps:
        if s not in STEPS:
            raise ValueError(f'Invalid skip step: {s!r} (valid steps: {STEPS})')
    skip_steps = set(skip_steps)
    if start_step is not None:
        if start_step not in STEPS:
            raise ValueError(
                f'Invalid start step: {start_step!r} (valid steps: {STEPS})'
            )
        skip_steps.update(STEPS[: STEPS.index(start_step)])
    if end_step is not None:
        if end_step not in STEPS:
            raise ValueError(f'Invalid end step: {end_step!r} (valid steps: {STEPS})')
        skip_steps.update(STEPS[STEPS.index(end_step) + 1 :])

    if basic_navigation:
        skip_steps.add('animate')

    # Standardise file paths
    root_path = os.path.expandvars(os.path.expanduser(root_path))
    if background_path is not None:
        background_path = os.path.expandvars(os.path.expanduser(background_path))

    log(f'Running MIRI pipeline')
    log(f'Root path: {root_path!r}', time=False)
    if skip_steps:
        log(
            f'Skipping steps: {", ".join([repr(s) for s in STEPS if s in skip_steps])}',
            time=False,
        )
    else:
        log(f'Running all pipeline steps', time=False)
    log(f'Desaturate: {desaturate!r}', time=False)
    if groups_to_use:
        log(f'Groups to keep: {groups_to_use!r}', time=False)
    if basic_navigation:
        log(f'Basic navigation: {basic_navigation!r}', time=False)
    log(f'Number of dithers: {ndither!r}', time=False)
    print()

    # Run pipeline steps...

    # TODO desaturation stuff

    if 'stage1' not in skip_steps:
        log('Running reduction stage 1')
        paths_in = sorted(glob.glob(os.path.join(root_path, 'stage0', '*.fits')))
        output_dir = os.path.join(root_path, 'stage1')
        args_list = [(p, output_dir, stage1_kwargs or {}) for p in paths_in]
        check_path(output_dir)

        log(f'Processing {len(args_list)} files...')
        runmany(
            reduction_detector1_fn,
            args_list,
            desc='stage1',
            **reduction_parallel_kwargs,
        )
        log('Reduction stage 1 step complete\n')

    if 'stage2' not in skip_steps:
        log('Running reduction stage 2')
        paths_in = sorted(glob.glob(os.path.join(root_path, 'stage1', '*.fits')))
        output_dir = os.path.join(root_path, 'stage2')
        args_list = [(p, output_dir, stage2_kwargs or {}) for p in paths_in]
        check_path(output_dir)

        log(f'Processing {len(args_list)} files...')
        runmany(
            reduction_spec2_fn, args_list, desc='stage2', **reduction_parallel_kwargs
        )
        log('Reduction stage 2 step complete\n')


def reduction_detector1_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    path_in, output_dir, kwargs = args
    Detector1Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)


def reduction_spec2_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    path_in, output_dir, kwargs = args
    Spec2Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)


def reduction_spec3_fn(path_in: str) -> None:
    pass


def main(*args):
    pass
    # TODO argparse stuff


if __name__ == '__main__':
    main()
