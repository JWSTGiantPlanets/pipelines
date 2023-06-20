#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full JWST reduction pipeline for MIRI MRS data, including the standard reduction from 
stage0 to stage3, and custom pipeline steps for additional cleaning and data 
visualisation.

See STEP_DESCRIPTIONS below for a description of each step in this pipeline.


Setup
=====
If you are running a reduction on Leicester's ALICE HPC then the pipeline should work 
out of the box without any additional setup - see the example job submission script 
below. If you are running the pipeline on another machine, you will need to set up the
correct environment to run the pipeline:

The `reduction` step of the pipeline requires the following environment variables to be
set to load the and cache the appropriate CRDS reference files ::

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
`/data/saturn/SATURN-15N/stage0`), then run the following command on the command line ::

    python3 miri_pipeline.py /data/saturn/SATURN-15N

or from Python ::

    import miri_pipeline
    miri_pipeline.run_pipeline('/data/saturn/SATURN-15N')

This will run the full pipeline, and output data files appropriate directories (e.g. 
`/data/saturn/SATURN-15N/stage3`, `/data/saturn/SATURN-15N/plots` etc.).

By default, the full pipeline is effectively run twice, first with the redisdual 
defringe step disabled, then with it enabled. If you only want defringed or 
non-defringed data, you can customise this behaviour with the `defringe` argument or the
`--defringe` or `--no-defringe` flags.

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

The walltime and memory requirements have been tested for the giant planet observations
(i.e. 4 dithers, ~5 groups etc.). If your data has more dithers/groups/integrations 
etc., you may need to increase the walltime and decrease the number of nodes (ppn).

To use this script you will need to:
- replace `py310` in the `conda activate` line to the name of your conda environment
- replace the two references to `/data/.../SATURN-15N` with the path to your data :: 

    #!/bin/bash
    #
    #PBS -N MIRI_Pipeline
    #PBS -l walltime=18:00:00
    #PBS -l vmem=80gb
    #PBS -l nodes=1:ppn=8
    
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
    python3 /data/nemesis/jwst/scripts/oliver/pipelines/miri_pipeline.py /data/nemesis/jwst/MIRI_IFU/Saturn_2022nov13/SATURN-15N --parallel
    
    # Change permissions on modified files so that other users can use them
    chmod -R --quiet 777 /data/nemesis/jwst/MIRI_IFU/Saturn_2022nov13/SATURN-15N
    chmod -R --quiet 777 $CRDS_PATH
"""
STEP_DESCRIPTIONS = """
- `remove_groups`: Remove groups from the data (for use in desaturating the data) [optional].
- `reduction`: Run the standard JWST reduction pipeline (stage0 to stage3).
- `navigate`: Navigate reduced files.
- `desaturate`: Desaturate data using cubes with fewer groups [optional].
- `flat`: Correct flat effects using synthetic flat field cubes.
- `despike`: Clean cubes by removing extreme outlier pixels.
- `background`: Subtract background from cubes [optional].
- `plot`: Generate quick look summary plots of data.
- `animate`: Generate summary animations of data.
"""
CLI_EXAMPLES = """examples:

# Print a help message, including documentation for each argument
python3 miri_pipeline.py -h

# Run the full pipeline, with all steps enabled
python3 miri_pipeline.py /data/saturn/SATURN-15N

# Run the full pipeline, without any desaturation
python3 miri_pipeline.py /data/saturn/SATURN-15N --no-desaturate

# Run the full pipeline, with only the defringed data
python3 miri_pipeline.py /data/saturn/SATURN-15N --defringe

# Run the pipeline, including background subtraction
python3 miri_pipeline.py /data/saturn/SATURN-15N --background_path /data/saturn/SATURN-BACKGROUND

# Run the pipeline, but stop before creating any visualisations
python3 miri_pipeline.py /data/saturn/SATURN-15N --end_step despike

# Re-run the pipeline, skipping the initial reduction steps
python3 miri_pipeline.py /data/saturn/SATURN-15N --start_step desaturate

# Run the pipeline, passing custom arguments to different steps
python3 miri_pipeline.py /data/saturn/SATURN-RINGS --kwargs '{"reduction": {"stages": [2, 3]}, "animation": {"radius_factor": 2.5}}'
"""
import argparse
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
from parallel_tools import runmany
from tools import KeepMissingDict, all_combinations, log

# Central record of the filepaths for various steps of the reduction pipeline
PATH_FITS = os.path.join(
    '{stage}',
    'd{dither}{fringe}{nav}',
    'Level3_ch{channel}-{band}_s3d{nav}.fits',
)
PATH_PLOT = os.path.join(
    'plots',
    '{stage}',
    'd{dither}{fringe}',
    '{channel}{abc}-Level3_ch{channel}-{band}_s3d.png',
)
PATH_ANIMATION = os.path.join(
    'animation',
    '{stage}',
    'dither_comparison{fringe}',
    '{channel}{abc}_animation.mp4',
)
PATH_STAGE3_RAW = PATH_FITS.format_map(KeepMissingDict(stage='stage3', nav=''))
PATH_DATA = PATH_FITS.format_map(KeepMissingDict(nav='_nav'))
PATH_NAVIGATED = PATH_DATA.format_map(KeepMissingDict(stage='stage3'))
PATH_DESATURATED = PATH_DATA.format_map(KeepMissingDict(stage='stage3_desaturated'))
PATH_FLAT = PATH_DATA.format_map(KeepMissingDict(stage='stage4_flat'))
PATH_DESPIKED = PATH_DATA.format_map(KeepMissingDict(stage='stage5_despike'))
PATH_BACKGROUND = PATH_DATA.format_map(KeepMissingDict(stage='stage6_background'))

DATA_STAGES = [
    'stage3',
    'stage3_desaturated',
    'stage4_flat',
    'stage5_despike',
    'stage6_background',
]

# Get a useful default flat data path that works natively on Leicester's ALICE HPC
_saturn_flat_path = os.path.join(
    'MIRI_IFU/Saturn_2022nov13/flat_field/merged',
    'ch{channel}-{band}_constructed_flat{fringe}.fits',
)
_path_options = [
    os.path.join('/data/nemesis/jwst', _saturn_flat_path),
    os.path.join(os.path.dirname(__file__), '..', 'jwst_data', _saturn_flat_path),
    os.path.join(os.path.dirname(__file__), 'flat_field'),
]
DEFAULT_FLAT_DATA_PATH = _path_options[-1]
for _p in _path_options:
    if os.path.exists(os.path.dirname(_p)):
        DEFAULT_FLAT_DATA_PATH = os.path.normpath(_p)
        break

# Arguments to format the paths above
CHANNELS = ['1', '2', '3', '4']
BANDS = ['short', 'medium', 'long']
FRINGES = ['', '_fringe']
CHANNEL_LENGTH_ALIASES = {'short': 'A', 'medium': 'B', 'long': 'C'}

# Pipeline constants
Step: TypeAlias = Literal[
    'remove_groups',
    'reduce',
    'navigate',
    'desaturate',
    'flat',
    'despike',
    'background',
    'plot',
    'animate',
]
STEPS: list[Step] = [
    'remove_groups',
    'reduce',
    'navigate',
    'desaturate',
    'flat',
    'despike',
    'background',
    'plot',
    'animate',
]


def run_pipeline(
    root_path: str,
    *,
    defringe: bool | Literal['both'] = 'both',
    desaturate: bool = True,
    groups_to_use: list[int] | None = None,
    background_path: str | None = None,
    ndither: int = 4,
    parallel: float | bool = False,
    flat_data_path: str = DEFAULT_FLAT_DATA_PATH,
    basic_navigation: bool = False,
    skip_steps: list[Step] | set[Step] | None = None,
    start_step: Step | None = None,
    end_step: Step | None = None,
    reduction_kwargs: dict[str, Any] | None = None,
    navigation_kwargs: dict[str, Any] | None = None,
    desaturation_kwargs: dict[str, Any] | None = None,
    flat_kwargs: dict[str, Any] | None = None,
    despike_kwargs: dict[str, Any] | None = None,
    background_kwargs: dict[str, Any] | None = None,
    plot_kwargs: dict[str, Any] | None = None,
    animation_kwargs: dict[str, Any] | None = None,
) -> None:
    """
    Run the full MIRI MRS reduction pipeline, including the standard reduction from
    stage0 to stage3, and custom pipeline steps for additional cleaning and data
    visualisation.

    Examples ::

        from miri_pipeline import run_pipeline

        # Run the full pipeline, with all steps enabled
        run_pipeline('/data/saturn/SATURN-15N')

        # Run the full pipeline, without any desaturation
        run_pipeline('/data/saturn/SATURN-15N', desaturate=False)

        # Run the pipeline, including background subtraction
        run_pipeline(
            '/data/saturn/SATURN-15N',
            background_path='data/saturn/SATURN-BACKGROUND',
        )

        # Run the pipeline, but stop before creating any visualisations
        run_pipeline('/data/saturn/SATURN-15N', end_step='despike')

        # Re-run the pipeline, skipping the initial reduction steps
        run_pipeline('/data/saturn/SATURN-15N', start_step='desaturate')

        # Run the pipeline, passing custom arguments to different steps
        run_pipeline(
            '/data/saturn/SATURN-RINGS',
            reduction_kwargs={'stages': [2, 3]},
            animation_kwargs={'radius_factor': 2.5},
        )

    Args:
        root_path: Path to the root directory containing the data. This directory should
            contain a subdirectory with the stage 0 data (i.e. `root_path/stage0`).
            Additional directories will be created automatically for each step of the
            pipeline (e.g. `root_path/stage1`, `root_path/plots` etc.).
        defringe: Toggle defringing of the data. If True, defringe data will be saved
            and processed. If `'both'` (the default), the pipeline will be run twice,
            first without defringing, then again with defringing enabled. Defringed data
            will be saved in separate directories (e.g. `root_path/stage3/d1_fringe`),
            so both sets of data will be available for further analysis.
        desaturate: Toggle desaturation of the data. If True, the `reduction` step will
            be run for different numbers of groups, which are then combined in the
            `desaturate` step to produce a desaturated data cube. If False, the
            `reduction` step will be run only once, with the full number of groups, and
            the `desaturate` step will be skipped (i.e. going straight to the flat field
            step).
        groups_to_use: List of groups to reduce and use for desaturating the data. If
            this is `None` (the default), then all available groups will be used. If
            `desaturate` is False, then this argument will be ignored. The data with all
            groups will always be included.
        background_path: Path to directory containing background data. For example, if
            your `root_path` is `/data/saturn/SATURN-15N`, the `background_path` may
            be `/data/saturn/SATURN-BACKGROUND`. Note that background subtraction will
            require the background data to be already reduced. If no `background_path`
            is specified (the default), then no background subtraction will be
            performed.
        ndither: Number of dithers for each observation. Defaults to 4.
        parallel: Fraction of CPU cores to use when multiprocessing. Set to 0 to run
            serially, 1 to use all cores, 0.5 to use half of the cores, etc. If enabled,
            multiprocessing is used in the `reduction`, `despike`, `plot` and `animate`
            steps.
        flat_data_path: Optionally specify custom path to the flat field data. This path
            should contain `{channel}`, `{band}` and `{fringe}` placeholders, which will
            be replaced by appropriate values for each channel, band and defringe
            setting.
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

        reduction_kwargs, navigation_kwargs, desaturation_kwargs, flat_kwargs,
        despike_kwargs, spx_kwargs, plot_kwargs, animation_kwargs: Arguments are passed
            to the corresponding functions for each step of the pipeline. These can
            therefore be used to override the default values for each step. See the
            documentation for each function for details on the arguments.
    """
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
    flat_data_path = os.path.expandvars(os.path.expanduser(flat_data_path))

    log('Running MIRI pipeline')
    log(f'Root path: {root_path!r}', time=False)
    if skip_steps:
        log(
            f'Skipping steps: {", ".join([repr(s) for s in STEPS if s in skip_steps])}',
            time=False,
        )
    else:
        log('Running all pipeline steps', time=False)
    log(f'Defringe: {defringe!r}', time=False)
    log(f'Desaturate: {desaturate!r}', time=False)
    if groups_to_use:
        log(f'Groups to keep: {groups_to_use!r}', time=False)
    if basic_navigation:
        log(f'Basic navigation: {basic_navigation!r}', time=False)
    log(f'Number of dithers: {ndither!r}', time=False)
    print()

    if defringe == 'both':
        kwargs = dict(
            desaturate=desaturate,
            groups_to_use=groups_to_use,
            background_path=background_path,
            ndither=ndither,
            parallel=parallel,
            flat_data_path=flat_data_path,
            basic_navigation=basic_navigation,
            navigation_kwargs=navigation_kwargs,
            desaturation_kwargs=desaturation_kwargs,
            flat_kwargs=flat_kwargs,
            despike_kwargs=despike_kwargs,
            background_kwargs=background_kwargs,
            plot_kwargs=plot_kwargs,
            animation_kwargs=animation_kwargs,
        )

        # First, run the entire pipeline with defringe=False
        # (do defringe=False first as it is much faster)
        run_pipeline(
            root_path,
            defringe=False,
            skip_steps=skip_steps,
            reduction_kwargs=reduction_kwargs,
            **kwargs,
        )

        # Then run the pipeline from the defringe step onwards with defringe=True
        # (no need to run the remove groups and reduction stages 1&2 again)
        skip_steps |= {'remove_groups'}
        reduction_kwargs = reduction_kwargs or {}
        reduction_kwargs.setdefault('stages', [2.5, 3])
        run_pipeline(
            root_path,
            defringe=True,
            skip_steps=skip_steps,
            reduction_kwargs=reduction_kwargs,
            **kwargs,
        )

        log('MIRI pipeline completed for both defringe settings')
        print()
        return

    # Run pipeline steps...

    if desaturate and 'remove_groups' not in skip_steps:
        log('Removing groups from data...')
        files = sorted(glob.glob(os.path.join(root_path, 'stage0', '*.fits')))
        # Allow customisation of remove_groups() groups_to_use argument with kwargs
        kw = {**dict(groups_to_use=groups_to_use), **(desaturation_kwargs or {})}
        for _p in tqdm.tqdm(files, desc='Removing groups'):
            remove_groups.remove_groups_from_file(_p, **kw)
        log('Remove groups step complete\n')

    # Create list of paths to reduce (both 'main' dataset and reduced groups)
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

    if 'reduce' not in skip_steps:
        log('Reducing data...')
        for _p in group_root_paths:
            # Allow customisation of reduction parallel kw with reduction_kwargs
            kw = {**dict(parallel=parallel), **(reduction_kwargs or {})}
            reduce_jwst_miri.JWSTReduction(
                _p,
                defringe=defringe,
                _run_immediately=True,
                **kw,
            )
        log('Reduction step complete\n')

    if 'navigate' not in skip_steps:
        log('Navigating data...')
        paths_to_navigate = []
        for _p in group_root_paths:
            paths_to_navigate.extend(
                get_stage_paths(_p, PATH_STAGE3_RAW, defringe, ndither)
            )
        paths_to_navigate = [_p for _p in paths_to_navigate if os.path.exists(_p)]
        navigate_jwst_observations.load_kernels()
        navigate_jwst_observations.navigate_multiple(
            *paths_to_navigate, basic=basic_navigation, **navigation_kwargs or {}
        )
        log('Navigation step complete\n')

    if desaturate and 'desaturate' not in skip_steps:
        log('Desaturating data...')
        paths_in_list = [
            get_stage_paths(_p, PATH_NAVIGATED, defringe, ndither)
            for _p in group_root_paths
        ]
        paths_out = get_stage_paths(root_path, PATH_DESATURATED, defringe, ndither)
        args_list = []
        for paths_in, path_out in zip(zip(*paths_in_list), paths_out):
            if os.path.exists(paths_in[0]):
                args_list.append((paths_in, path_out))

        log(f'Processing {len(args_list)} files...')
        for paths_in, paths_out in tqdm.tqdm(args_list, desc='Desaturating'):
            desaturate_data.replace_saturated(
                paths_in, paths_out, **desaturation_kwargs or {}
            )
        log('Desaturation step complete\n')

    if 'flat' not in skip_steps:
        log('Applying flat fields...')
        pattern_in = PATH_DESATURATED if desaturate else PATH_NAVIGATED
        log(f'Input file pattern: {pattern_in!r}', time=False)
        log(f'Flat data path: {flat_data_path!r}', time=False)
        paths = [
            (p_in, p_out, p_flat)
            for p_in, p_out, p_flat in zip(
                get_stage_paths(root_path, pattern_in, defringe, ndither),
                get_stage_paths(root_path, PATH_FLAT, defringe, ndither),
                get_stage_paths(None, flat_data_path, defringe, ndither),
            )
            if os.path.exists(p_in)
        ]

        log(f'Processing {len(paths)} files...')
        for p_in, p_out, p_flat in tqdm.tqdm(paths, desc='Flat fielding'):
            flat_field.apply_flat(p_in, p_out, p_flat, **flat_kwargs or {})
        log('Flat fielding step complete\n')

    if 'despike' not in skip_steps:
        log('Despiking data...')
        paths = [
            (p_in, p_out)
            for p_in, p_out in zip(
                get_stage_paths(root_path, PATH_FLAT, defringe, ndither),
                get_stage_paths(root_path, PATH_DESPIKED, defringe, ndither),
            )
            if os.path.exists(p_in)
        ]

        log(f'Processing {len(paths)} files...')
        if parallel:
            despike_kwargs = despike_kwargs or {}
            despike_kwargs['progress_bar'] = False
        despike_args = [(p_in, p_out, despike_kwargs) for p_in, p_out in paths]
        runmany(despike_fn, despike_args, parallel_frac=parallel, desc='Despiking')
        log('Despiking step complete\n')

    if (background_path is not None) and ('background' not in skip_steps):
        log('Subtracting background...')
        log(f'Background path: {background_path!r}', time=False)
        paths = [
            (p_in, p_out, p_bg)
            for p_in, p_out, p_bg in zip(
                get_stage_paths(root_path, PATH_DESPIKED, defringe, ndither),
                get_stage_paths(root_path, PATH_BACKGROUND, defringe, ndither),
                get_stage_paths(
                    background_path,
                    PATH_DESPIKED.format_map(KeepMissingDict(dither='1')),
                    defringe,
                    ndither,
                ),
            )
            if os.path.exists(p_in) and os.path.exists(p_bg)
        ]
        log(f'Processing {len(paths)} files...')
        for p_in, p_out, p_bg in tqdm.tqdm(paths, desc='Background subtracting'):
            background_subtraction.subtract_background(
                p_in, p_out, p_bg, **flat_kwargs or {}
            )
        log('Background step complete\n')

    if 'plot' not in skip_steps:
        log('Generating quick look plots...')
        # Reverse stages so that we get plots of the 'final' version generated first
        for stage in reversed(DATA_STAGES):
            log(f'Generating plots for stage {stage}...')
            template_in = PATH_DATA.format_map(KeepMissingDict(stage=stage))
            template_out = PATH_PLOT.format_map(KeepMissingDict(stage=stage))
            paths = []
            for _p in group_root_paths:
                paths.extend(
                    [
                        (p_in, p_out)
                        for p_in, p_out in zip(
                            get_stage_paths(_p, template_in, defringe, ndither),
                            get_stage_paths(_p, template_out, defringe, ndither),
                        )
                        if os.path.exists(p_in)
                    ]
                )

            log(f'Processing {len(paths)} files...')
            plot_args = [(p_in, p_out, plot_kwargs) for p_in, p_out in paths]
            runmany(
                plot_fn, plot_args, parallel_frac=parallel, desc=f'Plotting {stage}'
            )

        log('Plotting step complete\n')

    if 'animate' not in skip_steps:
        log('Generating animations...')
        # Use dict for paths here as multiple input paths map to the same output
        paths_dict: dict[str, list[str]] = {}
        for stage in reversed(DATA_STAGES):
            template_in = PATH_DATA.format_map(KeepMissingDict(stage=stage))
            template_out = PATH_ANIMATION.format_map(KeepMissingDict(stage=stage))
            for p_in, p_out in zip(
                get_stage_paths(root_path, template_in, defringe, ndither),
                get_stage_paths(root_path, template_out, defringe, ndither),
            ):
                if os.path.exists(p_in):
                    paths_dict.setdefault(p_out, []).append(p_in)

        log(f'{len(paths_dict)} animations to process...')
        animation_kwargs = animation_kwargs or {}
        animation_kwargs.setdefault('print_output', False)
        if parallel:
            animation_kwargs.setdefault('progress_bar', False)
        animation_args = [
            (p_out, p_ins, animation_kwargs) for p_out, p_ins in paths_dict.items()
        ]
        runmany(
            animation_fn,
            animation_args,
            parallel_frac=parallel,
            desc='Animating',
        )
        log('Animation step complete\n')

    log('MIRI pipeline complete')
    log(f'Root path: {root_path!r}', time=False)
    log(f'Defringe: {defringe!r}', time=False)
    log(f'Desaturate: {desaturate!r}', time=False)
    if skip_steps:
        log(
            f'Skipped steps: {", ".join([repr(s) for s in STEPS if s in skip_steps])}',
            time=False,
        )
    print()


def despike_fn(args: tuple[str, str, dict[str, Any] | None]):
    p_in, p_out, despike_kwargs = args
    despike_data.despike_cube(p_in, p_out, **despike_kwargs or {})


def plot_fn(args: tuple[str, str, dict[str, Any] | None]):
    p_in, p_out, plot_kwargs = args
    jwst_summary_plots.make_summary_plot(p_in, p_out, **plot_kwargs or {})


def animation_fn(args: tuple[str, list[str], dict[str, Any] | None]):
    p_out, p_ins, animation_kwargs = args
    jwst_summary_animation.make_animation(p_ins, p_out, **animation_kwargs or {})


def get_stage_paths(
    root_path: str | None,
    template: str,
    defringe: bool,
    ndither: int,
) -> list[str]:
    """
    Get a list of paths for all combinations of dithers, channels and bands for a given
    stage of the reduction.

    Args:
        root_path: Directory containing the data.
        template: Path of the data files relative to the `root_path`. The values
            `{dither}`, `{channel}`, `{band}`, `{abc}` and `{fringe}` will be replaced
            with the corresponding values for each individual file.
        defringe: Toggle defringing.
        ndither: Number of dithers.

    Returns:
        List of file paths, constructed from the root path and the template.
    """
    if root_path:
        template = os.path.join(root_path, template)
    fringe = '_fringe' if defringe else ''
    dithers = [str(i + 1) for i in range(ndither)]

    paths = []
    for path_kw in all_combinations(
        dither=dithers,
        channel=CHANNELS,
        band=BANDS,
    ):
        path_kw['abc'] = CHANNEL_LENGTH_ALIASES[path_kw['band']]
        paths.append(template.format(fringe=fringe, **path_kw))
    return paths


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description=(
            'Full JWST MIRI pipeline including the standard reduction from stage0 to '
            'stage3, and custom pipeline steps for additional cleaning and data '
            'visualisation. For more customisation, this script can be imported and '
            'run in Python (see the source code for mode details).\n\n'
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
        '--defringe',
        action=argparse.BooleanOptionalAction,
        default='both',
        help="""Toggle defringing of the data. If unspecified, the pipeline will be run 
            twice, first without defringing, then again with defringing enabled.""",
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
        '--groups_to_use',
        type=str,
        help="""Comma-separated list of groups to keep. For example, `1,2,3,4` will
            keep the first four groups. If unspecified, all groups will be kept.""",
    )
    parser.add_argument(
        '--background_path',
        type=str,
        help="""Path to directory containing background data. For example, if
            your `root_path` is `/data/saturn/SATURN-15N`, the `background_path` may
            be `/data/saturn/SATURN-BACKGROUND`. Note that background subtraction will
            require the background data to be already reduced. If no `background_path`
            is specified (the default), then no background subtraction will be 
            performed.""",
    )
    parser.add_argument(
        '--ndither',
        type=int,
        default=4,
        help="""Number of dithers in the dataset. Defaults to 4.""",
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
        equivalent to `--parallel 1`). If enabled, multiprocessing is used in the 
        `reduction`, `despike`, `plot` and `animate` steps.""",
    )
    parser.add_argument(
        '--flat_data_path',
        type=str,
        help="""Optionally specify custom path to the flat field data. This path
            should contain `{channel}`, `{band}` and `{fringe}` placeholders, which will
            be replaced by appropriate values for each channel, band and defringe
            setting.""",
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
            re-running part of the pipeline.""",
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
            `--kwargs '{"reduction": {"stages": [2, 3]}, "animation": {"radius_factor": 2.5}}'` 
            will pass the custom arguments to the reduction and animation steps.
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
