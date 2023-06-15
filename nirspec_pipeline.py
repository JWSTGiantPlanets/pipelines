#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NIRSPEC pipeline

TODO
- Remove groups (optional)
- Desaturation (optional)
- Add note about chmodding all the touched files
- Documentation
"""
import glob
import os
import pathlib
from typing import Any, Literal, TypeAlias

from astropy.io import fits
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline, Spec3Pipeline

import despike_data
import jwst_summary_animation
import jwst_summary_plots
import navigate_jwst_observations
from parallel_tools import runmany
from tools import check_path, log

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

DATA_STAGES = [
    'stage3',
    'stage3_desaturated',
    'stage4_despike',
]

# Default arguments for each step. These are merged with the user-supplied kwargs
# before being passed to the pipeline.
DEFAULT_STAGE1_KWARGS = {}
DEFAULT_STAGE2_KWARGS = {
    'steps': {'cube_build': {'coord_system': 'ifualign'}, 'extract_1d': {'skip': True}}
}
DEFAULT_STAGE3_KWARGS = {
    'steps': {'cube_build': {'coord_system': 'ifualign'}, 'extract_1d': {'skip': True}}
}
DEFAULT_NAVIGATION_KWARGS = {}
DEFAULT_DESPIKE_KWARGS = {}
DEFAULT_PLOT_KWARGS = {}
DEFAULT_ANIMATE_KWARGS = {}


def run_pipeline(
    root_path: str,
    *,
    desaturate: bool = True,
    groups_to_use: list[int] | None = None,
    parallel: float | bool = False,
    basic_navigation: bool = False,
    skip_steps: list[STEP] | set[STEP] | None = None,
    start_step: STEP | None = None,
    end_step: STEP | None = None,
    stage1_kwargs: dict[str, Any] | None = None,
    stage2_kwargs: dict[str, Any] | None = None,
    stage3_kwargs: dict[str, Any] | None = None,
    navigation_kwargs: dict[str, Any] | None = None,
    despike_kwargs: dict[str, Any] | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    animate_kwargs: dict[str, Any] | None = None,
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

    log('Running MIRI pipeline')
    log(f'Root path: {root_path!r}', time=False)
    if skip_steps:
        log(
            f'Skipping steps: {", ".join([repr(s) for s in STEPS if s in skip_steps])}',
            time=False,
        )
    else:
        log('Running all pipeline steps', time=False)
    log(f'Desaturate: {desaturate!r}', time=False)
    if groups_to_use:
        log(f'Groups to keep: {groups_to_use!r}', time=False)
    if basic_navigation:
        log(f'Basic navigation: {basic_navigation!r}', time=False)
    print()

    # Run pipeline steps...

    # TODO desaturation: remove groups
    # TODO desaturation: loop through stages1-3 with removed groups

    if 'stage1' not in skip_steps:
        run_stage1(root_path, stage1_kwargs, reduction_parallel_kwargs)

    if 'stage2' not in skip_steps:
        run_stage2(root_path, stage2_kwargs, reduction_parallel_kwargs)

    if 'stage3' not in skip_steps:
        run_stage3(root_path, stage3_kwargs, reduction_parallel_kwargs)

    if 'navigate' not in skip_steps:
        run_navigation(root_path, basic_navigation, navigation_kwargs)

    # TODO desaturation: combine removed groups data

    if 'despike' not in skip_steps:
        run_despike(root_path, despike_kwargs, parallel_kwargs)

    if 'plot' not in skip_steps:
        run_plot(root_path, plot_kwargs, parallel_kwargs)

    if 'animate' not in skip_steps:
        run_animate(root_path, animate_kwargs, parallel_kwargs)

    log('NIRSPEC pipeline complete')
    log(f'Root path: {root_path!r}', time=False)
    log(f'Desaturate: {desaturate!r}', time=False)
    if skip_steps:
        log(
            f'Skipped steps: {", ".join([repr(s) for s in STEPS if s in skip_steps])}',
            time=False,
        )
    print()


# stage1
def run_stage1(
    root_path: str,
    stage1_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
):
    log('Running reduction stage 1')
    paths_in = sorted(glob.glob(os.path.join(root_path, 'stage0', '*.fits')))
    output_dir = os.path.join(root_path, 'stage1')
    kwargs = DEFAULT_STAGE1_KWARGS | (stage1_kwargs or {})
    args_list = [(p, output_dir, kwargs) for p in paths_in]
    log(f'Output directory: {output_dir!r}', time=False)
    log(f'Arguments: {kwargs!r}', time=False)
    log(f'Processing {len(args_list)} files...', time=False)
    check_path(output_dir)
    runmany(
        reduction_detector1_fn,
        args_list,
        desc='stage1',
        **reduction_parallel_kwargs or {},
    )
    log('Reduction stage 1 step complete\n')


def reduction_detector1_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    path_in, output_dir, kwargs = args
    Detector1Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)


# stage2
def run_stage2(
    root_path: str,
    stage2_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running reduction stage 2')
    paths_in = sorted(glob.glob(os.path.join(root_path, 'stage1', '*rate.fits')))
    output_dir = os.path.join(root_path, 'stage2')
    kwargs = DEFAULT_STAGE2_KWARGS | (stage2_kwargs or {})
    args_list = [(p, output_dir, kwargs) for p in paths_in]
    log(f'Output directory: {output_dir!r}', time=False)
    log(f'Arguments: {kwargs!r}', time=False)
    log(f'Processing {len(args_list)} files...', time=False)
    check_path(output_dir)
    runmany(
        reduction_spec2_fn, args_list, desc='stage2', **reduction_parallel_kwargs or {}
    )
    log('Reduction stage 2 step complete\n')


def reduction_spec2_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    path_in, output_dir, kwargs = args
    Spec2Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)


# stage3
def run_stage3(
    root_path: str,
    stage3_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running reduction stage 3')
    kwargs = DEFAULT_STAGE3_KWARGS | (stage3_kwargs or {})
    log(f'Arguments: {kwargs!r}', time=False)

    grouped_files = group_stage2_files_for_stage3(root_path)
    for dither, paths_grouped in grouped_files.items():
        log(
            f'Processing dither {dither}/{len(grouped_files)} ({len(paths_grouped)} files)...'
        )
        output_dir = os.path.join(root_path, 'stage3', f'd{dither}')
        log(f'Output directory: {output_dir!r}', time=False)
        check_path(output_dir)
        asn_paths = []
        for filter_grating, paths in paths_grouped.items():
            asnfile = os.path.join(output_dir, f'l3asn-{filter_grating}.json')
            write_asn_for_stage3(paths, asnfile, prodname='Level3')
            asn_paths.append(asnfile)

        args_list = [(p, output_dir, stage3_kwargs or {}) for p in sorted(asn_paths)]

        runmany(
            reduction_spec3_fn,
            args_list,
            desc='stage3',
            **reduction_parallel_kwargs or {},
        )
    log('Reduction stage 3 step complete\n')


def group_stage2_files_for_stage3(root_path: str) -> dict[int, dict[str, list[str]]]:
    out: dict[int, dict[str, list[str]]] = {}
    paths_in = sorted(glob.glob(os.path.join(root_path, 'stage2', '*_cal.fits')))
    for p in paths_in:
        with fits.open(p) as hdul:
            hdr = hdul[0].header  #  type: ignore
        dither = int(hdr['PATT_NUM'])
        filter_ = hdr['FILTER']
        grating = hdr['GRATING']
        k = f'{filter_}_{grating}'
        out.setdefault(dither, {}).setdefault(k, []).append(p)
    return out


def write_asn_for_stage3(files: list, asnfile: str, prodname: str, **kwargs) -> None:
    asn = asn_from_list(files, rule=DMS_Level3_Base, product_name=prodname)
    if 'bg' in kwargs:
        for bgfile in kwargs['bg']:
            asn['products'][0]['members'].append(
                {'expname': bgfile, 'exptype': 'background'}
            )
    _, serialized = asn.dump()
    with open(asnfile, 'w', encoding='utf-8') as outfile:
        outfile.write(serialized)


def reduction_spec3_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    asn_path, output_dir, kwargs = args
    Spec3Pipeline.call(
        asn_path,
        output_dir=output_dir,
        save_results=True,
        **kwargs,
    )


# navigation
def run_navigation(
    root_path: str,
    basic_navigation: bool = False,
    navigation_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running navigation')
    kwargs = DEFAULT_NAVIGATION_KWARGS | (navigation_kwargs or {})
    paths_in = sorted(glob.glob(os.path.join(root_path, 'stage3', 'd*', '*_s3d.fits')))
    log(f'Basic navigation: {basic_navigation!r}', time=False)
    log(f'Arguments: {kwargs!r}', time=False)
    navigate_jwst_observations.load_kernels()
    navigate_jwst_observations.navigate_multiple(
        *paths_in, basic=basic_navigation, **kwargs
    )
    log('Navigation step complete\n')


# despike
def run_despike(
    root_path: str,
    despike_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
    input_stage: str = 'stage3',
) -> None:
    log('Running despike')
    paths_in = sorted(
        glob.glob(os.path.join(root_path, input_stage, 'd*_nav', '*_nav.fits'))
    )
    paths_out = [
        replace_path_part(p, -3, 'stage4_despike', check_old=input_stage)
        for p in paths_in
    ]
    parallel_kwargs = parallel_kwargs or {}
    kwargs = DEFAULT_DESPIKE_KWARGS | (despike_kwargs or {})
    if parallel_kwargs.get('parallel_frac', False):
        kwargs['progress_bar'] = False
    args_list = [(p_in, p_out, kwargs) for p_in, p_out in zip(paths_in, paths_out)]
    log(f'Arguments: {kwargs!r}', time=False)
    log(f'Processing {len(args_list)} files...', time=False)
    runmany(despike_fn, args_list, desc='despike', **parallel_kwargs)
    log('Despike step complete\n')


def despike_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    p_in, p_out, kwargs = args
    despike_data.despike_cube(p_in, p_out, **kwargs)


# plot
def run_plot(
    root_path: str,
    plot_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Generating quick look plots')
    kwargs = DEFAULT_PLOT_KWARGS | (plot_kwargs or {})
    parallel_kwargs = parallel_kwargs or {}
    log(f'Arguments: {kwargs!r}', time=False)
    for stage in reversed(DATA_STAGES):
        log(f'Generating plots for stage {stage!r}')
        paths_in = sorted(
            glob.glob(os.path.join(root_path, stage, 'd*_nav', '*_nav.fits'))
        )
        paths_out = [
            replace_path_part(p, -3, 'plots', check_old=stage) for p in paths_in
        ]
        paths_out = [
            replace_path_suffix(p, '.png', check_old='.fits') for p in paths_out
        ]
        args_list = [(p_in, p_out, kwargs) for p_in, p_out in zip(paths_in, paths_out)]
        log(f'Processing {len(args_list)} files...', time=False)
        runmany(plot_fn, args_list, desc='plot', **parallel_kwargs)
    log('Plot step complete\n')


def plot_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    p_in, p_out, kwargs = args
    jwst_summary_plots.make_summary_plot(p_in, p_out, **kwargs)


# animate
def run_animate(
    root_path: str,
    animate_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Generating quick look animations')
    kwargs = DEFAULT_ANIMATE_KWARGS | (animate_kwargs or {})
    parallel_kwargs = parallel_kwargs or {}
    log(f'Arguments: {kwargs!r}', time=False)
    for stage in reversed(DATA_STAGES):
        log(f'Generating animations for stage {stage!r}')
        paths_in = sorted(
            glob.glob(os.path.join(root_path, stage, 'd*_nav', '*_nav.fits'))
        )
        paths_dict: dict[str, list[str]] = {}
        for p_in in paths_in:
            filename = os.path.basename(
                replace_path_suffix(p_in, '.mp4', check_old='.fits')
            )
            p_out = os.path.join(root_path, 'animation', stage, filename)
            paths_dict.setdefault(p_out, []).append(p_in)
        args_list = [
            (paths_in, p_out, kwargs) for p_out, paths_in in paths_dict.items()
        ]
        log(f'Processing {len(args_list)} files...', time=False)
        runmany(animation_fn, args_list, desc='animate', **parallel_kwargs)
    log('Animate step complete\n')


def animation_fn(args: tuple[list[str], str, dict[str, Any]]) -> None:
    paths_in, p_out, kwargs = args
    jwst_summary_animation.make_animation(paths_in, p_out, **kwargs)


# utils
def replace_path_part(
    path: str, idx: int, new: str, *, check_old: str | None = None
) -> str:
    """
    Replace a part of a path with a new value.

    Args:
        path: The path to modify.
        idx: The index of the part to replace in the directory tree.
        new: The new value to use.
        check_old: If not None, check that the value at the index is this value
            before replacing it.

    Returns:
        The modified path.
    """
    p = pathlib.Path(path)
    parts = list(p.parts)
    if check_old is not None and parts[idx] != check_old:
        raise ValueError(f'Expected {check_old!r} at index {idx} in {path!r}')
    parts[idx] = new
    return str(pathlib.Path(*parts))


def replace_path_suffix(path: str, new: str, *, check_old: str | None = None) -> str:
    """
    Replace the suffix of a path with a new value.

    Args:
        path: The path to modify.
        new: The new value to use.
        check_old: If not None, check that the suffix is this value before
            replacing it.

    Returns:
        The modified path.
    """
    p = pathlib.Path(path)
    if check_old is not None and p.suffix != check_old:
        raise ValueError(f'Expected {check_old!r} suffix in {path!r}')
    return str(p.with_suffix(new))


def main():
    pass
    # TODO argparse stuff


if __name__ == '__main__':
    main()
