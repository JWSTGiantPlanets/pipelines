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

Data cubes are saved at each step of the pipeline, with the data in each stage directory
used as the input for the next stage. The `stage0`, `stage1` and `stage2` directories
contain partially reduced data cubes. The `stage3` directory onwards contain reduced 
data cubes which can be used for science. Within each science directory (`stage3`, 
`stage4...` etc.), files are organised by dither and any data variant. For example:
- `stage3/d1` contains the stage3 data for dither 1
- `stage3/d2_bg` contains the stage3 data for dither 2 with the background subtracted
- `stage3/combined` contains the stage4 data with all dithers combined
- `stage4_despike/d3` contains the stage4 for dither 3

The `plots` and `animation` directories contain quick look visualisations of the data.

Metadata about the pipeline processing steps is saved in the `PRIMARY` FITS header of 
each file.

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
    chmod -R --quiet ugo+rw /data/uranus/lon1
    chmod -R --quiet ugo+rw $CRDS_PATH
"""
STEP_DESCRIPTIONS = """
- `remove_groups`: Remove groups from the data (for use in desaturating the data) [optional].
- `stage1`: Run the standard JWST reduction pipeline stage 1.
- `stage2`: Run the standard JWST reduction pipeline stage 2 (including optional background subtraction).
- `stage3`: Run the standard JWST reduction pipeline stage 3.
- `navigate`: Navigate reduced files.
- `desaturate`: Desaturate data using cubes with fewer groups [optional].
- `despike`: Clean cubes by removing extreme outlier pixels.
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
python3 nirspec_pipeline.py /data/uranus/lon1 --kwargs '{"stage3": {"steps": {"outlier_detection": {"snr": "30.0 24.0", "scale": "1.3 0.7"}}}, "plot": {"plot_brightest_spectrum": true}}'
"""
from typing import Any, Collection

from jwst.assign_wcs.util import NoDataOnDetectorError

from pipeline import BoolOrBoth, Pipeline, Step, get_pipeline_argument_parser

# Pipeline constants for NIRSpec
STEPS: tuple[Step, ...] = (
    'remove_groups',
    'stage1',
    'stage2',
    'stage3',
    'navigate',
    'desaturate',
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
        }
    },
    'stage3': {
        'steps': {
            'extract_1d': {'skip': True},
            'cube_build': {'coord_system': 'ifualign'},
        }
    },
}
STEP_DIRECTORIES: dict[Step, tuple[str, str]] = {'despike': ('', 'stage4_despike')}
STAGE_DIRECTORIES_TO_PLOT = (
    'stage3',
    'stage3_desaturated',
    'stage4_despike',
)


def run_pipeline(
    root_path: str,
    *,
    skip_steps: Collection[Step] | None = None,
    start_step: Step | None = None,
    end_step: Step | None = None,
    parallel: float | bool = False,
    desaturate: bool = True,
    groups_to_use: list[int] | None = None,
    background_subtract: BoolOrBoth = 'both',
    background_path: str | None = None,
    basic_navigation: bool = False,
    step_kwargs: dict[Step, dict[str, Any]] | None = None,
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
    """
    pipeline = NirspecPipeline(
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
    )
    pipeline.run(
        skip_steps=skip_steps,
        start_step=start_step,
        end_step=end_step,
    )


class NirspecPipeline(Pipeline):
    @staticmethod
    def get_instrument() -> str:
        return 'NIRSpec'

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
        directories['despike'] = (dir_in, directories['despike'][1])
        return directories

    @property
    def stage3_file_match_hdr_keys(self) -> tuple[str, ...]:
        return ('FILTER', 'GRATING')

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        return STAGE_DIRECTORIES_TO_PLOT

    # Step overrides
    def reduction_spec2_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        try:
            super().reduction_spec2_fn(args)
        except NoDataOnDetectorError:
            path_in, output_dir, kwargs = args
            print(f'No data on detector for {path_in!r}, skipping')


def main():
    parser = get_pipeline_argument_parser(
        NirspecPipeline, STEP_DESCRIPTIONS, CLI_EXAMPLES
    )
    run_pipeline(**vars(parser.parse_args()))


if __name__ == '__main__':
    main()
