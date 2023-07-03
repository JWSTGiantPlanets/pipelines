#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full JWST reduction pipeline for MIRI MRS data, including the standard reduction from 
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
`/data/uranus/lon1/stage0`), then run the following command on the command line ::

    python3 miri_pipeline.py /data/uranus/lon1

or from Python ::

    import miri_pipeline
    miri_pipeline.run_pipeline('/data/uranus/lon1')

This will run the full pipeline, and output data files appropriate directories (e.g. 
`/data/uranus/lon1/stage3`, `/data/uranus/lon1/plots` etc.).

By default, most steps of the pipeline are effectively run twice, each step with the 
residual defringe step disabled first, then with it enabled. If you only want defringed
or non-defringed data, you can customise this behaviour with the `defringe` argument or 
the `--defringe` or `--no-defringe` flags.

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

The walltime and memory requirements have been tested for the giant planet observations
(i.e. 4 dithers, ~5 groups etc.). If your data has more dithers/groups/integrations 
etc., you may need to increase the walltime and decrease the number of nodes (ppn).

To use this script you will need to:
- replace `py310` in the `conda activate` line to the name of your conda environment
- replace the two references to `/data/uranus/lon1` with the path to your data :: 

    #!/bin/bash
    #
    #PBS -N Pipeline
    #PBS -l walltime=24:00:00
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
    python3 /data/nemesis/jwst/scripts/oliver/pipelines/miri_pipeline.py /data/uranus/lon1 --parallel
    
    # Change permissions on modified files so that other users can use them
    chmod -R --quiet 777 /data/uranus/lon1
    chmod -R --quiet 777 $CRDS_PATH
"""
STEP_DESCRIPTIONS = """
- `remove_groups`: Remove groups from the data (for use in desaturating the data) [optional].
- `stage1`: Run the standard JWST reduction pipeline stage 1.
- `stage2`: Run the standard JWST reduction pipeline stage 2 (including optional background subtraction).
- `defringe`: Run the JWST reduction pipeline residual fringe step [optional].
- `stage3`: Run the standard JWST reduction pipeline stage 3.
- `navigate`: Navigate reduced files.
- `desaturate`: Desaturate data using cubes with fewer groups [optional].
- `flat`: Correct flat effects using synthetic flat field cubes.
- `despike`: Clean cubes by removing extreme outlier pixels.
- `plot`: Generate quick look summary plots of data.
- `animate`: Generate summary animations of data.
"""
CLI_EXAMPLES = """examples:

# Print a help message, including documentation for each argument
python3 miri_pipeline.py -h

# Run the full pipeline, with all steps enabled
python3 miri_pipeline.py /data/uranus/lon1

# Run the full pipeline, without any desaturation
python3 miri_pipeline.py /data/uranus/lon1 --no-desaturate

# Run the full pipeline, with only the defringed data
python3 miri_pipeline.py /data/uranus/lon1 --defringe

# Run the pipeline, including background subtraction
python3 miri_pipeline.py /data/uranus/lon1 --background_path /data/uranus/background

# Run the pipeline, but stop before creating any visualisations
python3 miri_pipeline.py /data/uranus/lon1 --end_step despike

# Only run the plotting step
python3 miri_pipeline.py /data/uranus/lon1 --start_step plot --end_step plot

# Re-run the pipeline, skipping the initial reduction steps
python3 miri_pipeline.py /data/uranus/lon1 --start_step desaturate

# Run the pipeline, passing custom arguments to different steps
python3 miri_pipeline.py /data/uranus/lon1 --kwargs '{"stage3": {"steps": {"outlier_detection": {"snr": "30.0 24.0", "scale": "1.3 0.7"}}}, "plot": {"plot_brightest_spectrum": true}}'
"""
import argparse
import os
from typing import Any, Collection, Literal

import tqdm
from astropy.io import fits
from jwst.residual_fringe import ResidualFringeStep

import flat_field
from parallel_tools import runmany
from pipeline import Pipeline, Step, get_pipeline_argument_parser

# Pipeline constants for MIRI
STEPS = (
    'remove_groups',
    'stage1',
    'stage2',
    'defringe',
    'stage3',
    'navigate',
    'desaturate',
    'flat',
    'despike',
    'plot',
    'animate',
)
DEFAULT_KWARGS: dict[Step, dict[str, Any]] = {
    'stage2': {
        'steps': {
            'cube_build': {'skip': True},
            'extract_1d': {'skip': True},
            'bkg_subtract': {'skip': False},
        },
    },
    'stage3': {
        'steps': {
            'extract_1d': {'skip': True},
            'cube_build': {'output_type': 'band', 'coord_system': 'ifualign'},
        }
    },
}
STAGE_DIRECTORIES_TO_PLOT = (
    'stage3',
    'stage3_desaturated',
    'stage4_flat',
    'stage5_despike',
    'stage6_background',
)
STEP_DIRECTORIES: dict[Step, tuple[str, str]] = {
    'defringe': ('stage2', 'stage2'),
    'flat': ('', 'stage4_flat'),
    'despike': ('stage4_flat', 'stage5_despike'),
}

BAND_ABC_ALIASES = {'short': 'A', 'medium': 'B', 'long': 'C'}


# Try to get a useful default flat data path that at least works natively on Leicester's
# ALICE HPC
_path_options = (
    '/data/nemesis/jwst/MIRI_IFU/flat_field',
    os.path.join(
        os.path.dirname(__file__), '..', 'jwst_data', 'MIRI_IFU', 'flat_field'
    ),
    os.path.join(os.path.dirname(__file__), '..', 'jwst_data', 'flat_field'),
    os.path.join(os.path.dirname(__file__), 'flat_field'),
)
DEFAULT_FLAT_DATA_PATH = _path_options[-1]
for _p in _path_options:
    if os.path.exists(os.path.dirname(_p)):
        DEFAULT_FLAT_DATA_PATH = os.path.normpath(_p)
        break


def run_pipeline(
    root_path: str,
    *,
    skip_steps: Collection[Step] | None = None,
    start_step: Step | None = None,
    end_step: Step | None = None,
    parallel: float | bool = False,
    desaturate: bool = True,
    groups_to_use: list[int] | None = None,
    background_subtract: bool | Literal['both'] = 'both',
    background_path: str | None = None,
    basic_navigation: bool = False,
    step_kwargs: dict[Step, dict[str, Any]] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
    defringe: bool | Literal['both'] = 'both',
    flat_data_path: str = DEFAULT_FLAT_DATA_PATH,
) -> None:
    """
    Run the full MIRI IFU reduction pipeline, including the standard reduction from
    stage0 to stage3, and custom pipeline steps for additional cleaning and data
    visualisation.

    Examples ::

        from miri_pipeline import run_pipeline

        # Run the full pipeline with all steps enabled
        run_pipeline('/data/uranus/lon1')

        # Run the full pipeline in parallel, using 50% of available cores
        run_pipeline('/data/uranus/lon1', parallel=0.5)

        # Run the full pipeline, without any desaturation
        run_pipeline('/data/uranus/lon1', desaturate=False)

        # Run the pipeline, including background subtraction in stage2
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
            kwargs{
                'stage3: {
                    'outlier_detection': {'snr': '30.0 24.0', 'scale': '1.3 0.7'}
                },
                'plot': {'plot_brightest_spectrum': True},
                'animate': {'radius_factor': 2.5},
            },
        )

    Args:
        root_path: Path to the root directory containing the data. This directory should
            contain a subdirectory with the stage 0 data (i.e. `root_path/stage0`).
            Additional directories will be created automatically for each step of the
            pipeline (e.g. `root_path/stage1`, `root_path/plots` etc.). This is the only
            required argument.

        skip_steps: List of steps to skip.
        start_step: Convenience argument to add all steps before `start_step` to
            `skip_steps`.
        end_step: Convenience argument to add all steps after `end_step` to
            `skip_steps`.

        parallel: Fraction of CPU cores to use when multiprocessing. Set to 0 to run
            serially, 1 to use all cores, 0.5 to use half of the cores, etc.
        desaturate: Toggle desaturation of the data. If True, the `stage1`-`navigate`
            steps will be run for different numbers of groups, which are then combined
            in the`desaturate` step to produce a desaturated data cube.
        groups_to_use: List of groups to reduce and use for desaturating the data. If
            this is `None` (the default), then all available groups will be used. If
            `desaturate` is False, then this argument will be ignored. The data with all
            groups will always be included.
        background_subtract: Toggle background subtraction. If True, the backgrounds
            in the `background_path` directory will be subtracted from the data in
            `stage2`. If 'both' (the default), then versions with and without background
            subtraction will be created.
        background_path: Path to the directory containing the background data (i.e. the
            equivalent of `root_path` for the background observations). Note that the
            background data must be reduced to at least `stage1` for this to work as
            the *_rate.fits files are used for the background subtraction.
        basic_navigation: Toggle between basic or full navigation. If True, then only
            RA and Dec navigation backplanes are generated (e.g. useful for small
            bodies). If False (the default), then full navigation is performed,
            generating a full set of coordinate backplanes (lon/lat, illumination
            angles etc.). Using basic navigation automatically skips the animation step.
            This is mainly useful if you get SPICE errors when navigating the data.
        step_kwargs: Dictionary of keyword arguments to pass to each step of the
            pipeline. The keys should be the step names (e.g. `stage1`, `stage2` etc.)
            and the values should be dictionaries of keyword arguments to pass to the
            step. See the documentation for each step for the available keyword
            arguments. For example,
            `step_kwargs={'stage3':{'outlier_detection': {'snr': '30.0 24.0', 'scale': '1.3 0.7'}}}`
            can will customise the `stage3` step. See the DEFAULT_KWARGS constant above
            for the default values.
        parallel_kwargs: Dictionary of keyword arguments to customise parallel
            processing.
        reduction_parallel_kwargs: Dictionary of keyword arguments to customise parallel
            processing for the reduction steps (i.e. `stage1`-`stage3`). This will be
            merged with `parallel_kwargs` (i.e.
            `parallel_kwargs | reduction_parallel_kwargs`).

        defringe: Toggle defringing of the data. If True, defringe data will be saved
            and processed. If `'both'` (the default), the pipeline steps will
            effectively be run twice, first without defringing, then again with
            defringing enabled. Defringed data will be saved in separate directories
            (e.g. `root_path/stage3/d1_fringe`), so both sets of data will be available
            for further analysis.
        flat_data_path: Optionally specify custom path to the flat field data. This path
            should contain `{channel}`, `{band}` and `{fringe}` placeholders, which will
            be replaced by appropriate values for each channel, band and defringe
            setting.
    """
    pipeline = MiriPipeline(
        root_path,
        parallel=parallel,
        desaturate=desaturate,
        groups_to_use=groups_to_use,
        background_subtract=background_subtract,
        background_path=background_path,
        basic_navigation=basic_navigation,
        step_kwargs=step_kwargs,
        parallel_kwargs=parallel_kwargs,
        reduction_parallel_kwargs=reduction_parallel_kwargs,
        flat_data_path=flat_data_path,
        defringe=defringe,
    )
    pipeline.run(
        skip_steps=skip_steps,
        start_step=start_step,
        end_step=end_step,
    )


class MiriPipeline(Pipeline):
    def __init__(
        self,
        *args,
        flat_data_path: str,
        defringe: bool | Literal['both'] = 'both',
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.flat_data_path = self.standardise_path(flat_data_path)
        self.defringe = defringe

    @staticmethod
    def get_instrument() -> str:
        return 'MIRI'

    @property
    def steps(self) -> tuple[Step, ...]:
        return STEPS

    @property
    def default_kwargs(self) -> dict[Step, dict[str, Any]]:
        return DEFAULT_KWARGS

    @property
    def step_directories(self) -> dict[Step, tuple[str, str]]:
        directories = super().step_directories | STEP_DIRECTORIES

        dir_in = 'stage3_desaturated' if self.desaturate else 'stage3'
        directories['flat'] = (dir_in, directories['flat'][1])

        return directories

    @property
    def stage3_file_match_hdr_keys(self) -> tuple[str, ...]:
        return ('CHANNEL', 'BAND')

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        return STAGE_DIRECTORIES_TO_PLOT

    def process_steps_list(
        self,
        skip_steps: Collection[Step] | None,
        start_step: Step | None,
        end_step: Step | None,
    ) -> set[Step]:
        skip_steps = super().process_steps_list(skip_steps, start_step, end_step)
        if not self.defringe:
            skip_steps.add('defringe')
        return skip_steps

    @property
    def data_variants_individual(self) -> set[str]:
        variants = super().data_variants_individual
        if self.defringe:
            variants.add('fringe')
        return variants

    def print_reduction_info(self, skip_steps: set[Step]) -> None:
        super().print_reduction_info(skip_steps)
        self.log(f'Defringe: {self.defringe!r}', time=False)
        self.log(f'Flat data path: {self.flat_data_path!r}', time=False)

    # Step overrides
    def reduction_detector1_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        kwargs = kwargs.copy()
        with fits.open(path_in) as hdul:
            ngroups = hdul['PRIMARY'].header['NGROUPS']  # type: ignore
        if ngroups <= 3:
            kwargs['steps'] = kwargs.get('steps', {}) | {
                'firstframe': {'skip': True},
                'lastframe': {'skip': True},
                'rscd': {'skip': True},
                'jump': {'skip': True},
            }
        return super().reduction_detector1_fn((path_in, output_dir, kwargs))

    def run_defringe(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['defringe']
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(root_path, dir_in, '*cal.fits')
            if 'bg' in self.data_variants_individual:
                paths_in += self.get_paths(root_path, dir_in, 'bg', '*cal.fits')
            args_list = [(p, os.path.basename(p), kwargs) for p in paths_in]
            self.log(f'Processing {len(args_list)} files...', time=False)
            runmany(
                self.defringe_fn,
                args_list,
                desc='defringe',
                **self.reduction_parallel_kwargs,
            )

    def defringe_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        ResidualFringeStep.call(
            path_in, output_dir=output_dir, save_results=True, **kwargs
        )

    def get_stage3_variant_paths_in(self, variant: frozenset[str]) -> list[str]:
        """
        Get list of input paths for a given variant for stage3.
        """
        dir_in, dir_out = self.step_directories['stage3']
        if variant == frozenset():
            return self.get_paths(dir_in, '*cal.fits')
        elif variant == frozenset({'bg'}):
            return self.get_paths(dir_in, 'bg', '*cal.fits')
        elif variant == frozenset({'bg', 'fringe'}):
            return self.get_paths(dir_in, 'bg', '*residual_fringe.fits')
        elif variant == frozenset({'fringe'}):
            return self.get_paths(dir_in, '*residual_fringe.fits')
        raise ValueError(f'Unknown variant: {variant}')

    def run_flat(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['flat']
        self.log(f'Flat data path: {self.flat_data_path!r}', time=False)
        # Use d* to skip dither combined files
        paths_in = self.get_paths(dir_in, 'd*', '*_nav.fits', filter_variants=True)
        self.log(f'Processing {len(paths_in)} files...', time=False)
        for p_in in tqdm.tqdm(paths_in, desc='flat'):
            p_out = self.replace_path_part(p_in, -3, dir_out, old=dir_in)
            with fits.open(p_in) as hdul:
                hdr = hdul['PRIMARY'].header  #  type: ignore
            p_flat = self.flat_data_path.format(
                channel=hdr['CHANNEL'],
                band=hdr['BAND'].lower(),
                fringe='_fringe' if hdr['S_RESFRI'] == 'COMPLETE' else '',
            )
            flat_field.apply_flat(p_in, p_out, p_flat, **kwargs)

    def get_plot_filename_prefix(self, path: str) -> str:
        with fits.open(path) as hdul:
            hdr = hdul['PRIMARY'].header  #  type: ignore
        channel = hdr['CHANNEL']
        abc = BAND_ABC_ALIASES[hdr['BAND'].casefold().strip()]
        return f'{channel}{abc}'


def main():
    parser = get_pipeline_argument_parser(MiriPipeline, STEP_DESCRIPTIONS, CLI_EXAMPLES)
    parser.add_argument(
        '--defringe',
        action=argparse.BooleanOptionalAction,
        default='both',
        help="""Toggle defringing of the data. If unspecified, the pipeline steps will
            effectively be run twice, first without defringing, then again with
            defringing enabled.""",
    )
    parser.add_argument(
        '--flat_data_path',
        type=str,
        help="""Optionally specify custom path to the flat field data. This path
            should contain `{channel}`, `{band}` and `{fringe}` placeholders, which will
            be replaced by appropriate values for each channel, band and defringe
            setting.""",
    )
    run_pipeline(**vars(parser.parse_args()))


if __name__ == '__main__':
    main()
