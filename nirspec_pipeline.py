#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full JWST reduction pipeline for NIRSpec IFU data, including the standard reduction from
stage0 to stage3, and custom pipeline steps for additional cleaning and data 
visualisation.

See STEP_DESCRIPTIONS below for a description of each step in this pipeline.

GitHub repository: https://github.com/JWSTGiantPlanets/pipelines


Setup
=====
If you are running a reduction on Leicester's ALICE HPC then the pipeline should work 
out of the box without any additional setup - see the example job submission script 
below. If you are running the pipeline on another machine, you will need to set up the
correct environment to run the pipeline:

The `stage1`-`stage3` steps of the pipeline require the following environment variables
to be set to load the and cache the appropriate CRDS reference files ::

    export CRDS_PATH="path/to/crds_cache"
    export CRDS_SERVER_URL="https://jwst-crds.stsci.edu"

If you are starting with an empty CRDS cache, the pipeline will have to download several
GB of reference files, so the first run will be much slower. If you are getting errors
when the reduction pipeline is reading/writing CRDS reference files, try running the 
reduction step serially (i.e. without the `--parallel` flag) to avoid any potential race
conditions caused by multiple processes trying to simultaneously download the same 
reference files.

The navigation and visualisation steps use SPICE kernels to calculate the observation
geometry. The pipeline will attempt to automatically find the SPICE kernels, but this 
will generally only work if they are saved in `~/spice_kernels` or you are running the
code on ALICE. Therefore, the location of these kernels can customised by setting the
environment variable ::

    export PLANETMAPPER_KERNEL_PATH="path/to/spice_kernels"

For more information on downloading and saving the SPICE kernels, see 
https://planetmapper.readthedocs.io/en/latest/spice_kernels.html.


Usage
=====
The pipeline can be run from the command line, or imported and run from Python. To run
the pipeline to fully reduce a dataset, simply download the stage0 (e.g. to 
`/data/uranus/lon1/stage0`), then run the following command on the command line ::

    python3 nirspec_pipeline.py /data/uranus/lon1

or from Python ::

    import nirspec_pipeline
    nirspec_pipeline.run_pipeline('/data/uranus/lon1')

This will run the full pipeline, and output data files appropriate directories (e.g. 
`/data/uranus/lon1/stage3`, `/data/uranus/lon1/plots` etc.).

If the pipeline is run with desaturation enabled, the pipeline flow is:
- Firstly, multiple versions of the stage0 cubes are created with different numbers of 
  groups (the `remove_groups` step).
- `stage1` is run on e.g. the 4 group data, then the 3 group data, then the 2 group
  data, then the 1 group data. Then `stage2` is run on the 4-1 group data, then 
  `stage3`, then `navigate`.
- The `desaturate` step is then run to combine the various group data into a single
  desaturated dataset.
- Subsequent pipeline steps (e.g. `plot`) are run on the desaturated data.

For more command line examples, see CLI_EXAMPLES below, and for more Python examples,
see the docstring for `run_pipeline()` below.


Customising logging
===================
The the `reduction` step of the pipeline can print a large amount of information to the 
terminal. To prevent this, you can create a `stpipe-log.cfg` file in your working 
directory, with the following contents to redirect the output to a file named 
`pipeline.log` (this `pipeline.log` file will be created automatically if needed, so you
don't need to create it yourself) ::

    [*]
    handler = append:pipeline.log
    level = WARNING

See https://jwst-pipeline.readthedocs.io/en/latest/jwst/stpipe/user_logging.html for
more details.


Example ALICE HPC job submission script
=======================================
The following example job submission script will run the pipeline on ALICE, using
multiprocessing for the `reduction` step. The chmod commands at the end change the 
permissions on data and CRDS cache files so that any new/modified files can be accessed 
by other users in the future.

Depending on the number of dithers/groups/integrations etc. for your data,
you may need to increase the walltime and decrease the number of nodes (ppn).

To use this script you will need to:
- replace `py310` in the `conda activate` line to the name of your conda environment
- replace the two references to `/data/uranus/lon1` with the path to your data :: 

    #!/bin/bash
    #
    #PBS -N NIRSpec_Pipeline
    #PBS -l walltime=24:00:00
    #PBS -l vmem=80gb
    #PBS -l nodes=1:ppn=4
    
    source ~/.bashrc 
    conda activate py310 # <-- replace this with the name of your conda environment

    export CRDS_PATH="/data/nemesis/jwst/crds_cache"
    export CRDS_SERVER_URL="https://jwst-crds.stsci.edu"

    # Optionally redirect the verbose `reduction` step output to the file `pipeline.log`
    if [ -f "/data/nemesis/jwst/scripts/oliver/pipelines/stpipe-log.cfg" ]; then
        cp -f "/data/nemesis/jwst/scripts/oliver/pipelines/stpipe-log.cfg" .
        echo "Copied stpipe-log.cfg to current working directory"
    fi

    # Run the pipeline
    python3 /data/nemesis/jwst/scripts/oliver/pipelines/nirspec_pipeline.py /data/uranus/lon1 --parallel
    
    # Change permissions on modified files so that other users can use them
    chmod -R --quiet 777 /data/uranus/lon1
    chmod -R --quiet 777 $CRDS_PATH
"""
STEP_DESCRIPTIONS = """
- `remove_groups`: Remove groups from the data (for use in desaturating the data) [optional].
- `stage1`: Run the standard JWST reduction pipeline stage 1.
- `stage2`: Run the standard JWST reduction pipeline stage 2.
- `stage3`: Run the standard JWST reduction pipeline stage 3.
- `navigate`: Navigate reduced files.
- `desaturate`: Desaturate data using cubes with fewer groups [optional].
- `despike`: Clean cubes by removing extreme outlier pixels.
- `background`: Subtract background from cubes [optional].
- `plot`: Generate quick look summary plots of data.
- `animate`: Generate summary animations of data.
"""
CLI_EXAMPLES = """examples:

# Print a help message, including documentation for each argument
python3 nirspec_pipeline.py -h

# Run the full pipeline with all steps enabled
python3 nirspec_pipeline.py /data/uranus/lon1

# Run the full pipeline in parallel using all cores
python3 nirspec_pipeline.py /data/uranus/lon1 --parallel

# Run the full pipeline in parallel, using 50% of available cores
python3 nirspec_pipeline.py /data/uranus/lon1 --parallel 0.5

# Run the full pipeline, without any desaturation
python3 nirspec_pipeline.py /data/uranus/lon1 --no-desaturate

# Run the pipeline, including background subtraction
python3 nirspec_pipeline.py /data/uranus/lon1 --background_path /data/uranus/background

# Run the pipeline, but stop before creating any visualisations
python3 nirspec_pipeline.py /data/uranus/lon1 --end_step despike

# Only run the plotting step
python3 nirspec_pipeline.py /data/uranus/lon1 --start_step plot --end_step plot

# Re-run the pipeline, skipping the initial reduction steps
python3 nirspec_pipeline.py /data/uranus/lon1 --start_step desaturate

# Run the pipeline, passing custom arguments to different steps
python3 nirspec_pipeline.py /data/uranus/lon1 --kwargs '{"stage3": {"steps": {"outlier_detection": {"snr": "30.0 24.0", "scale": "1.3 0.7"}}}, "plot": {"plot_brightest_spectrum": true}, "animation": {"radius_factor": 2.5}}'
"""
import argparse
import glob
import json
import os
import pathlib
from typing import Any, Literal, TypeAlias

import tqdm
from astropy.io import fits
from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline, Spec3Pipeline

import background_subtraction
import desaturate_data
import despike_data
import jwst_summary_animation
import jwst_summary_plots
import navigate_jwst_observations
import remove_groups
from parallel_tools import runmany
from tools import check_path, log

Step: TypeAlias = Literal[
    'remove_groups',
    'stage1',
    'stage2',
    'stage3',
    'navigate',
    'desaturate',
    'despike',
    'background',
    'plot',
    'animate',
]
STEPS: list[Step] = [
    'remove_groups',
    'stage1',
    'stage2',
    'stage3',
    'navigate',
    'desaturate',
    'despike',
    'background',
    'plot',
    'animate',
]

DATA_STAGES = [
    'stage3',
    'stage3_desaturated',
    'stage4_despike',
    'stage5_background',
]

# Default arguments for each step. These are merged with the user-supplied kwargs
# before being passed to the pipeline, for example:
# stage1_kwargs = DEFAULT_STAGE1_KWARGS | stage1_kwargs
DEFAULT_STAGE1_KWARGS = {}
DEFAULT_STAGE2_KWARGS = {
    'steps': {'cube_build': {'coord_system': 'ifualign'}, 'extract_1d': {'skip': True}}
}
DEFAULT_STAGE3_KWARGS = {
    'steps': {'cube_build': {'coord_system': 'ifualign'}, 'extract_1d': {'skip': True}}
}
DEFAULT_DESATURATE_KWARGS = {}
DEFAULT_NAVIGATE_KWARGS = {}
DEFAULT_DESPIKE_KWARGS = {}
DEFAULT_BACKGROUND_KWARGS = {}
DEFAULT_PLOT_KWARGS = {}
DEFAULT_ANIMATE_KWARGS = {}


def run_pipeline(
    root_path: str,
    *,
    desaturate: bool = True,
    groups_to_use: list[int] | None = None,
    background_path: str | None = None,
    parallel: float | bool = False,
    basic_navigation: bool = False,
    skip_steps: list[Step] | set[Step] | None = None,
    start_step: Step | None = None,
    end_step: Step | None = None,
    stage1_kwargs: dict[str, Any] | None = None,
    stage2_kwargs: dict[str, Any] | None = None,
    stage3_kwargs: dict[str, Any] | None = None,
    desaturate_kwargs: dict[str, Any] | None = None,
    navigate_kwargs: dict[str, Any] | None = None,
    despike_kwargs: dict[str, Any] | None = None,
    background_kwargs: dict[str, Any] | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    animate_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    """
    Run the full NIRSpec IFU reduction pipeline, including the standard reduction from
    stage0 to stage3, and custom pipeline steps for additional cleaning and data
    visualisation.

    Examples ::

        from nirspec_pipeline import run_pipeline

        # Run the full pipeline with all steps enabled
        run_pipeline('/data/uranus/lon1')

        # Run the full pipeline in parallel, using 50% of available cores
        run_pipeline('/data/uranus/lon1', parallel=0.5)

        # Run the full pipeline, without any desaturation
        run_pipeline('/data/uranus/lon1', desaturate=False)

        # Run the pipeline, including background subtraction
        run_pipeline(
            '/data/uranus/lon1',
            background_path='data/uranus/background',
        )

        # Run the pipeline, but stop before creating any visualisations
        run_pipeline('/data/uranus/lon1', end_step='despike')

        # Re-run the pipeline, skipping the initial reduction steps
        run_pipeline('/data/uranus/lon1', start_step='desaturate')

        # Run the pipeline, passing custom arguments to different steps
        run_pipeline(
            '/data/uranus/lon1',
            stage3_kwargs={
                'outlier_detection': {'snr': '30.0 24.0', 'scale': '1.3 0.7'}
            },
            plot_kwargs={'plot_brightest_spectrum': True},
            animation_kwargs={'radius_factor': 2.5},
        )

    Args:
        root_path: Path to the root directory containing the data. This directory should
            contain a subdirectory with the stage 0 data (i.e. `root_path/stage0`).
            Additional directories will be created automatically for each step of the
            pipeline (e.g. `root_path/stage1`, `root_path/plots` etc.).
        desaturate: Toggle desaturation of the data. If True, the `stage1` to `navigate`
            steps will be run for different numbers of groups, which are then combined
            in the `desaturate` step to produce a desaturated data cube. If False, the
            `reduction` step will be run only once, with the full number of groups, and
            the `desaturate` step will be skipped (i.e. going straight to the despike
            step).
        parallel: Fraction of CPU cores to use when multiprocessing. Set to 0 to run
            serially, 1 to use all cores, 0.5 to use half of the cores, etc.
        groups_to_use: List of groups to reduce and use for desaturating the data. If
            this is `None` (the default), then all available groups will be used. If
            `desaturate` is False, then this argument will be ignored. The data with all
            groups will always be included.
        background_path: Path to directory containing background data. For example, if
            your `root_path` is `/data/uranus/lon1`, the `background_path` may
            be `/data/uranus/background`. Note that background subtraction will
            require the background data to be already reduced. If no `background_path`
            is specified (the default), then no background subtraction will be
            performed.
        basic_navigation: Toggle between basic or full navigation. If True, then only
            RA and Dec navigation backplanes are generated (e.g. useful for small
            bodies). If False (the default), then full navigation is performed,
            generating a full set of coordinate backplanes (lon/lat, illumination
            angles etc.). Using basic navigation automatically skips the animation step.
        skip_steps: List of steps to skip. Defaults to None. See `STEPS` above for a
            list of valid steps. This is mainly useful if you are re-running part of the
            pipeline.
        start_step: Convenience argument to add all steps before `start_step` to
            `skip_steps`.
        end_step: Convenience argument to add all steps after `end_step` to
            `skip_steps`.

        stage1_kwargs, stage2_kwargs, stage3_kwargs, desaturate_kwargs,
        navigate_kwargs, despike_kwargs, background_kwargs, plot_kwargs, animate_kwargs:
            Dictionaries of arguments passed to the corresponding functions for each
            step of the pipeline. These can therefore be used to override the default
            values for each step. See the documentation for each function for details on
            the arguments. See the constants (e.g. DEFAULT_STAGE1_KWARGS) above for the
            default values.
    """
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

    # Process skip steps
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

    # Print info
    log('Running NIRSpec pipeline')
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
    if desaturate and 'remove_groups' not in skip_steps:
        run_remove_groups(root_path, groups_to_use)

    # Create list of paths to reduce (both 'main' dataset and reduced groups)
    group_root_paths = get_group_root_paths(root_path, desaturate, groups_to_use)

    if 'stage1' not in skip_steps:
        run_stage1(group_root_paths, stage1_kwargs, reduction_parallel_kwargs)

    if 'stage2' not in skip_steps:
        run_stage2(group_root_paths, stage2_kwargs, reduction_parallel_kwargs)

    if 'stage3' not in skip_steps:
        run_stage3(group_root_paths, stage3_kwargs, reduction_parallel_kwargs)

    if 'navigate' not in skip_steps:
        run_navigate(group_root_paths, basic_navigation, navigate_kwargs)

    if desaturate and 'desaturate' not in skip_steps:
        run_desaturate(root_path, group_root_paths, desaturate_kwargs, parallel_kwargs)

    if 'despike' not in skip_steps:
        run_despike(
            root_path,
            despike_kwargs,
            parallel_kwargs,
            input_stage='stage3_desaturated' if desaturate else 'stage3',
        )

    if background_path is not None and 'background' not in skip_steps:
        run_background(root_path, background_path, background_kwargs)

    if 'plot' not in skip_steps:
        run_plot(group_root_paths, plot_kwargs, parallel_kwargs)

    if 'animate' not in skip_steps:
        run_animate(root_path, animate_kwargs, parallel_kwargs)

    log('NIRSpec pipeline completed with')
    log(f'Root path: {root_path!r}', time=False)
    log(f'Desaturate: {desaturate!r}', time=False)
    if skip_steps:
        log(
            f'Skipped steps: {", ".join([repr(s) for s in STEPS if s in skip_steps])}',
            time=False,
        )
    log('ALL DONE')


# remove groups
def run_remove_groups(root_path: str, groups_to_use: list[int] | None = None) -> None:
    log('Removing groups')
    paths_in = sorted(glob.glob(os.path.join(root_path, 'stage0', '*_uncal.fits')))
    log(f'Processing {len(paths_in)} files...', time=False)
    for p in tqdm.tqdm(paths_in, desc='remove_groups'):
        remove_groups.remove_groups_from_file(p, groups_to_use)
    log('Remove groups step complete\n')


def get_group_root_paths(
    root_path: str, desaturate: bool, groups_to_use: list[int] | None = None
) -> list[str]:
    group_root_paths = [root_path]
    if desaturate:
        # If desaturating, we also need to reduce the data with fewer groups
        # This list is sorted such that the number of groups is decreasing (so that the
        # desaturation works correctly)
        reduced_group_root_paths = sorted(
            glob.glob(os.path.join(root_path, 'groups', '*_groups')),
            reverse=True,
            key=lambda _p: int(os.path.basename(_p).split('_')[0]),
        )
        if groups_to_use is not None:
            reduced_group_root_paths = [
                _p
                for _p in reduced_group_root_paths
                if int(os.path.basename(_p).split('_')[0]) in groups_to_use
            ]
        group_root_paths.extend(reduced_group_root_paths)
    return group_root_paths


# stage1
def run_stage1(
    group_root_paths: list[str],
    stage1_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
):
    log('Running reduction stage 1')
    kwargs = DEFAULT_STAGE1_KWARGS | (stage1_kwargs or {})
    log(f'Arguments: {kwargs!r}', time=False)
    for root_path in group_root_paths:
        if len(group_root_paths) > 1:
            log(f'Running reduction stage 1 for {root_path!r}')
        paths_in = sorted(glob.glob(os.path.join(root_path, 'stage0', '*.fits')))
        output_dir = os.path.join(root_path, 'stage1')
        args_list = [(p, output_dir, kwargs) for p in paths_in]
        log(f'Output directory: {output_dir!r}', time=False)
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
    group_root_paths: list[str],
    stage2_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running reduction stage 2')
    kwargs = DEFAULT_STAGE2_KWARGS | (stage2_kwargs or {})
    log(f'Arguments: {kwargs!r}', time=False)
    for root_path in group_root_paths:
        if len(group_root_paths) > 1:
            log(f'Running reduction stage 2 for {root_path!r}')
        paths_in = sorted(glob.glob(os.path.join(root_path, 'stage1', '*rate.fits')))
        output_dir = os.path.join(root_path, 'stage2')
        args_list = [(p, output_dir, kwargs) for p in paths_in]
        log(f'Output directory: {output_dir!r}', time=False)
        log(f'Processing {len(args_list)} files...', time=False)
        check_path(output_dir)
        runmany(
            reduction_spec2_fn,
            args_list,
            desc='stage2',
            **reduction_parallel_kwargs or {},
        )
    log('Reduction stage 2 step complete\n')


def reduction_spec2_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    path_in, output_dir, kwargs = args
    try:
        Spec2Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)
    except NoDataOnDetectorError:
        log(f'No data on detector for {path_in!r}, skipping')


# stage3
def run_stage3(
    group_root_paths: list[str],
    stage3_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running reduction stage 3')
    kwargs = DEFAULT_STAGE3_KWARGS | (stage3_kwargs or {})
    log(f'Arguments: {kwargs!r}', time=False)

    for root_path in group_root_paths:
        if len(group_root_paths) > 1:
            log(f'Running reduction stage 3 for {root_path!r}')

        grouped_files = group_stage2_files_for_stage3(root_path)

        # Only need to include the tile in prodname if it is needed to avoid filename
        # collisions for observations which use mosaicing. Most observations just have a
        # single tile per datset, so we can just use the standard prodname in this case
        keys_with_tiles = set(grouped_files.keys())
        keys_without_tiles = set(
            (dither, filter_, grating)
            for dither, tile, filter_, grating in keys_with_tiles
        )
        include_tile_in_prodname = len(keys_with_tiles) != len(keys_without_tiles)

        asn_paths_list: list[tuple[str, str]] = []
        for (dither, tile, filter_, grating), paths in grouped_files.items():
            dirname = 'combined' if dither is None else f'd{dither}'
            output_dir = os.path.join(root_path, 'stage3', dirname)
            check_path(output_dir)
            asn_path = os.path.join(
                output_dir, f'l3asn-{tile}_{filter_}_{grating}_dither-{dirname}.json'
            )
            write_asn_for_stage3(
                paths,
                asn_path,
                prodname='Level3' + (f'_{tile}' if include_tile_in_prodname else ''),
            )
            asn_paths_list.append((asn_path, output_dir))
        args_list = [
            (p, output_dir, kwargs) for p, output_dir in sorted(asn_paths_list)
        ]
        log(f'Processing {len(args_list)} files...', time=False)
        runmany(
            reduction_spec3_fn,
            args_list,
            desc='stage3',
            **reduction_parallel_kwargs or {},
        )
    log('Reduction stage 3 step complete\n')


def group_stage2_files_for_stage3(
    root_path: str,
) -> dict[tuple[int | None, str, str, str], list[str]]:
    out: dict[tuple[int | None, str, str, str], list[str]] = {}
    paths_in = sorted(glob.glob(os.path.join(root_path, 'stage2', '*_cal.fits')))
    for p in paths_in:
        with fits.open(p) as hdul:
            hdr = hdul[0].header  #  type: ignore
        dither = int(hdr['PATT_NUM'])
        tile = hdr['ACT_ID']
        filter_ = hdr['FILTER']
        grating = hdr['GRATING']

        # Keep dithers separate
        k = (dither, tile, filter_, grating)
        out.setdefault(k, []).append(p)

        # Combine dithers
        k = (None, tile, filter_, grating)
        out.setdefault(k, []).append(p)
    return out


def write_asn_for_stage3(
    files: list[str], asnfile: str, prodname: str, **kwargs
) -> None:
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
def run_navigate(
    group_root_paths: list[str],
    basic_navigation: bool = False,
    navigate_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running navigate')
    kwargs = DEFAULT_NAVIGATE_KWARGS | (navigate_kwargs or {})
    log(f'Basic navigation: {basic_navigation!r}', time=False)
    log(f'Arguments: {kwargs!r}', time=False)
    navigate_jwst_observations.load_kernels()
    for root_path in group_root_paths:
        if len(group_root_paths) > 1:
            log(f'Running navigation for {root_path!r}')
        paths_in = sorted(
            glob.glob(os.path.join(root_path, 'stage3', '*', '*_s3d.fits'))
        )
        navigate_jwst_observations.navigate_multiple(
            *paths_in, basic=basic_navigation, **kwargs
        )
    log('Navigate step complete\n')


# desaturate
def run_desaturate(
    root_path: str,
    group_root_paths: list[str],
    desaturate_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running desaturate')
    kwargs = DEFAULT_DESATURATE_KWARGS | (desaturate_kwargs or {})
    parallel_kwargs = parallel_kwargs or {}
    log(f'Arguments: {kwargs!r}', time=False)

    paths_dict: dict[str, list[str]] = {}
    for rp in group_root_paths:
        paths_in = sorted(glob.glob(os.path.join(rp, 'stage3', '*_nav', '*_nav.fits')))
        for p_in in paths_in:
            relpath = os.path.relpath(p_in, rp)
            relpath = replace_path_part(
                relpath, -3, 'stage3_desaturated', check_old='stage3'
            )
            p_out = os.path.join(root_path, relpath)
            paths_dict.setdefault(p_out, []).append(p_in)
    args_list = [(paths_in, p_out, kwargs) for p_out, paths_in in paths_dict.items()]
    log(f'Processing {len(args_list)} output files...', time=False)
    runmany(desaturate_fn, args_list, desc='desaturate', **parallel_kwargs)
    log('Desaturate step complete\n')


def desaturate_fn(args: tuple[list[str], str, dict[str, Any]]) -> None:
    paths_in, path_out, kwargs = args
    desaturate_data.replace_saturated(paths_in, path_out, **kwargs)


# despike
def run_despike(
    root_path: str,
    despike_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
    input_stage: str = 'stage3',
) -> None:
    log('Running despike')
    paths_in = sorted(
        glob.glob(os.path.join(root_path, input_stage, '*_nav', '*_nav.fits'))
    )
    paths_out = [
        replace_path_part(p, -3, 'stage4_despike', check_old=input_stage)
        for p in paths_in
    ]
    parallel_kwargs = parallel_kwargs or {}
    kwargs = DEFAULT_DESPIKE_KWARGS | (despike_kwargs or {})
    if parallel_kwargs.get('parallel_frac', False):
        kwargs.setdefault('progress_bar', False)
    args_list = [(p_in, p_out, kwargs) for p_in, p_out in zip(paths_in, paths_out)]
    log(f'Arguments: {kwargs!r}', time=False)
    log(f'Processing {len(args_list)} files...', time=False)
    runmany(despike_fn, args_list, desc='despike', **parallel_kwargs)
    log('Despike step complete\n')


def despike_fn(args: tuple[str, str, dict[str, Any]]) -> None:
    p_in, p_out, kwargs = args
    despike_data.despike_cube(p_in, p_out, **kwargs)


# background
def run_background(
    root_path: str,
    background_path: str,
    background_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running background')
    kwargs = DEFAULT_BACKGROUND_KWARGS | (background_kwargs or {})
    log(f'Arguments: {kwargs!r}', time=False)

    paths_in = sorted(
        glob.glob(os.path.join(root_path, 'stage4_despike', 'd*_nav', '*_nav.fits'))
    )
    log(f'Processing {len(paths_in)} input files...', time=False)
    for p_in in tqdm.tqdm(paths_in, desc='background'):
        p_out = replace_path_part(
            p_in, -3, 'stage5_background', check_old='stage4_despike'
        )
        with fits.open(p_in) as hdul:
            hdr = hdul['PRIMARY'].header  #  type: ignore
        p_background = os.path.join(
            background_path,
            'stage4_despike',
            'd{dither}_nav',
            'Level3_{grating}-{filter}_s3d_nav.fits',
        ).format(
            grating=hdr['GRATING'].lower(),
            filter=hdr['FILTER'].lower(),
            dither='1',
        )
        background_subtraction.subtract_background(p_in, p_out, p_background, **kwargs)
    log('Background step complete\n')


# plot
def run_plot(
    group_root_paths: list[str],
    plot_kwargs: dict[str, Any] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
) -> None:
    log('Running plot')
    kwargs = DEFAULT_PLOT_KWARGS | (plot_kwargs or {})
    parallel_kwargs = parallel_kwargs or {}
    log(f'Arguments: {kwargs!r}', time=False)
    for root_path in group_root_paths:
        if len(group_root_paths) > 1:
            log(f'Generating quick look plots for {root_path!r}')
        for stage in reversed(DATA_STAGES):
            log(f'Generating quick look plots for stage {stage!r}')
            paths_in = sorted(
                glob.glob(os.path.join(root_path, stage, '*_nav', '*_nav.fits'))
            )
            paths_out = [
                replace_path_part(p, -3, os.path.join('plots', stage), check_old=stage)
                for p in paths_in
            ]
            paths_out = [
                replace_path_suffix(p, '.png', check_old='.fits') for p in paths_out
            ]
            args_list = [
                (p_in, p_out, kwargs) for p_in, p_out in zip(paths_in, paths_out)
            ]
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
    log('Running animate')
    kwargs = DEFAULT_ANIMATE_KWARGS | (animate_kwargs or {})
    parallel_kwargs = parallel_kwargs or {}
    if parallel_kwargs.get('parallel_frac', False):
        kwargs.setdefault('progress_bar', False)
    log(f'Arguments: {kwargs!r}', time=False)
    for stage in reversed(DATA_STAGES):
        log(f'Generating quick look animations for stage {stage!r}')
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
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            'Full JWST NIRSpec IFU pipeline including the standard reduction from '
            'stage0 to stage3, and custom pipeline steps for additional cleaning and '
            'data visualisation. For more customisation, this script can be imported '
            'and run in Python (see the source code for mode details).\n\n'
            'The following steps are run in the full pipeline:' + STEP_DESCRIPTIONS
        ),
        epilog=CLI_EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        argument_default=argparse.SUPPRESS,
    )
    parser.add_argument(
        'root_path',
        type=str,
        help="""Path to the directory containing the data. This directory should
            contain a subdirectory with the stage0 data (i.e. `root_path/stage0`).
            Additional directories will be created automatically for each step of the
            pipeline (e.g. `root_path/stage1`, `root_path/plots` etc.).""",
    )
    parser.add_argument(
        '--desaturate',
        action=argparse.BooleanOptionalAction,
        default=True,
        help="""Toggle desaturation of the data. If desaturation is enabled the 
            `reduction` step will be run for different numbers of groups, which are then
            combined in the `desaturate` step to produce a desaturated data cube. This
            desaturation is enabled by default.
            """,
    )
    parser.add_argument(
        '--parallel',
        nargs='?',
        const=1,
        type=float,
        help="""Fraction of CPU cores to use when multiprocessing. For example, set 
        0.5 to use half of the available cores. If unspecified, then multiprocessing 
        will not be used and the pipeline will be run serially. If specified but no 
        value is given, then all available cores will be used (i.e. `--parallel` is 
        equivalent to `--parallel 1`).""",
    )
    parser.add_argument(
        '--groups_to_use',
        type=str,
        help="""Comma-separated list of groups to keep. For example, `1,2,3,4` will
            keep the first four groups. If unspecified, all groups will be kept.""",
    )
    parser.add_argument(
        '--background_path',
        type=str,
        help="""Path to directory containing background data. For example, if
            your `root_path` is `/data/uranus/lon1`, the `background_path` may
            be `/data/uranus/background`. Note that background subtraction will
            require the background data to be already reduced. If no `background_path`
            is specified (the default), then no background subtraction will be 
            performed.""",
    )
    parser.add_argument(
        '--basic_navigation',
        action='store_true',
        help="""Use basic navigation, and only save RA and Dec backplanes (e.g. useful
            for small bodies). By default, full navigation is performed, generating a
            full set of coordinate backplanes (lon/lat, illumination angles etc.). Using
            basic navigation automatically skips the navigation step.""",
    )
    parser.add_argument(
        '--skip_steps',
        nargs='+',
        type=str,
        help="""List of steps to skip. This is generally only useful if you are
            re-running part of the pipeline. Multiple steps can be passed as a
            space-separated list. For example, `--skip_steps flat despike`.""",
    )
    parser.add_argument(
        '--start_step',
        type=str,
        help="""Convenience argument to add all steps before `start_step` to 
            `skip_steps`.""",
    )
    parser.add_argument(
        '--end_step',
        type=str,
        help="""Convenience argument to add all steps steps after `end_step` to 
            `skip_steps`.""",
    )
    parser.add_argument(
        '--kwargs',
        type=str,
        help="""JSON string containing keyword arguments to pass to individual pipeline
            steps. For example, 
            `--kwargs '{"stage3": {"steps": {"outlier_detection": {"snr": "30.0 24.0", "scale": "1.3 0.7"}}}, "animation": {"radius_factor": 2.5}}'` 
            will pass the custom arguments to the stage3 and animation steps.
            """,
    )
    args = vars(parser.parse_args())
    json_kwargs = args.pop('kwargs', None)
    if json_kwargs:
        json_kwargs = json.loads(json_kwargs)
        for k, v in json_kwargs.items():
            if not k.endswith('_kwargs'):
                k = k + '_kwargs'
            args[k] = v
    if 'groups_to_use' in args:
        args['groups_to_use'] = [int(g) for g in args['groups_to_use'].split(',')]
    run_pipeline(**args)


if __name__ == '__main__':
    main()
