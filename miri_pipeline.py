#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full JWST reduction pipeline for MIRI MRS data, including the standard reduction from 
stage0 to stage3, and custom pipeline steps for additional cleaning and data 
visualisation.

See STEP_DESCRIPTIONS below for a description of each step in this pipeline.

GitHub repository: https://github.com/JWSTGiantPlanets/pipelines


Reference
=========
King et al., (2023). Custom JWST NIRSpec/IFU and MIRI/MRS Data Reduction Pipelines for
Solar System Targets. Research Notes of the AAS, 7(10), 223,
https://doi.org/10.3847/2515-5172/ad045f


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
or non-defringed data, you can customise this behaviour with the pipeline arguments/CLI
flags.

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
- `stage3/d2_bg_fringe` contains the stage3 data for dither 2 with the background
   subtracted and the residual fringes step run on the data
- `stage3/combined` contains the stage4 data with all dithers combined
- `stage4_flat/d3` contains the stage4 for dither 3

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

The walltime and memory requirements have been tested for the giant planet observations
(i.e. 4 dithers, ~5 groups etc.). If your data has more dithers/groups/integrations 
etc., you may need to increase the walltime and decrease the number of nodes (ppn).

To use this script you will need to:
- replace `py311` in the `source ~/miniconda...` line to your conda environemnt name
- replace the two references to `/data/uranus/lon1` with the path to your data :: 

    #!/bin/bash
    #
    #SBATCH --job-name=MIRI_Pipeline
    #SBATCH --time=24:00:00
    #SBATCH --mem=80G
    #SBATCH --ntasks=8
    #SBATCH --export=NONE
    #SBATCH --account=nemesis
    #SBATCH --output=slurm-%j-%x.out

    source ~/miniconda3/bin/activate py311

    export CRDS_PATH="/data/nemesis/jwst/crds_cache"
    export CRDS_SERVER_URL="https://jwst-crds.stsci.edu"
    export MPLBACKEND=Agg
    
    # Optionally redirect the verbose stage1-3 logging to the file `pipeline.log`
    if [ -f "/data/nemesis/jwst/scripts/oliver/pipelines/stpipe-log.cfg" ]; then
        cp -f "/data/nemesis/jwst/scripts/oliver/pipelines/stpipe-log.cfg" .
        echo "Copied stpipe-log.cfg to current working directory"
    fi

    # Run the pipeline
    python3 /data/nemesis/jwst/scripts/oliver/pipelines/miri_pipeline.py /data/uranus/lon1 --parallel
    exit_code=$?

    # Change permissions on modified files so that other users can use them
    chmod -R --quiet ugo+rw /data/uranus/lon1
    chmod -R --quiet ugo+rw $CRDS_PATH

    # Exit script with the same exit code as the pipeline
    exit $exit_code
"""
STEP_DESCRIPTIONS = """
- `remove_groups`: Remove groups from the data (for use in desaturating the data) [optional].
- `stage1`: Run the standard JWST reduction pipeline stage 1.
- `stage2`: Run the standard JWST reduction pipeline stage 2 (including optional background subtraction).
- `defringe`: Run the JWST reduction pipeline 2D residual fringe step (note `defringe_1d` is generally preferred) [optional].
- `stage3`: Run the standard JWST reduction pipeline stage 3 (including automatically running the spectral leak step on x1d spectra, if present).
- `defringe_1d`: Apply the JWST reduction pipeline 1D residual fringe step to each spectrum in the s3d data cubes [optional].
- `navigate`: Navigate reduced files.
- `psf`: Correct PSF effects using a forward model convolved with telescope's PSF [optional].
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

# Run the full pipeline, with desaturation
python3 miri_pipeline.py /data/uranus/lon1 --desaturate

# Run the full pipeline, with emsm weighting when building spectral cubes
python3 miri_pipeline.py /data/uranus/lon1 --cube-build-weighting=emsm

# Run the full pipeline, with only the defringed data
python3 miri_pipeline.py /data/uranus/lon1 --defringe_1d

# Run the pipeline, including background subtraction
python3 miri_pipeline.py /data/uranus/lon1 --background-path /data/uranus/background

# Run the pipeline, but stop before creating any visualisations
python3 miri_pipeline.py /data/uranus/lon1 --end-step despike

# Only run the plotting step
python3 miri_pipeline.py /data/uranus/lon1 --start-step plot --end-step plot

# Re-run the pipeline, skipping the initial reduction steps
python3 miri_pipeline.py /data/uranus/lon1 --desaturate --start-step desaturate

# Run the pipeline, passing custom arguments to different steps
python3 miri_pipeline.py /data/uranus/lon1 --kwargs '{"stage3": {"steps": {"outlier_detection": {"snr": "30.0 24.0", "scale": "1.3 0.7"}}}, "plot": {"plot_brightest_spectrum": true}}'
"""
import os
import pathlib
from typing import Any, Collection, Literal

import tqdm
from astropy.io import fits
from jwst.residual_fringe import ResidualFringeStep
from jwst.spectral_leak import SpectralLeakStep

import argparse_utils
import flat_field
import psf_correction
from defringe_1d import defringe_file
from parallel_tools import runmany
from pipeline import BoolOrBoth, Pipeline, RootPath, Step, get_pipeline_argument_parser

# Pipeline constants for MIRI
STEPS: tuple[Step, ...] = (
    'remove_groups',
    'stage1',
    'stage2',
    'defringe',
    'stage3',
    'defringe_1d',
    'navigate',
    'psf',
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
            'cube_build': {'output_type': 'band', 'coord_system': 'ifualign'},
        }
    },
}
STAGE_DIRECTORIES_TO_PLOT = (
    'stage3',
    'stage3_desaturated',
    'stage4_flat',
    'stage5_despike',
)
STEP_DIRECTORIES: dict[Step, tuple[str, str]] = {
    'defringe': ('stage2', 'stage2'),
    'defringe_1d': ('stage3', 'stage3'),
    'flat': ('', 'stage4_flat'),
    'psf': ('stage3', 'stage3'),
    'despike': ('stage4_flat', 'stage5_despike'),
}

BAND_ABC_ALIASES = {'short': 'A', 'medium': 'B', 'long': 'C'}


# Try to get a useful default flat data path that at least works natively on Leicester's
# ALICE HPC
DEFAULT_FLAT_FILENAME = 'ch{channel}-{band}_constructed_flat{fringe}.fits'
_path_options = (
    os.path.join('/data/nemesis/jwst/MIRI_IFU/flat_field', DEFAULT_FLAT_FILENAME),
    os.path.join(
        os.path.dirname(__file__),
        '..',
        'jwst_data',
        'MIRI_IFU',
        'flat_field',
        DEFAULT_FLAT_FILENAME,
    ),
    os.path.join(os.path.dirname(__file__), '..', 'jwst_data', 'flat_field', DEFAULT_FLAT_FILENAME),
    os.path.join(os.path.dirname(__file__), 'flat_field', DEFAULT_FLAT_FILENAME),
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
    desaturate: bool = False,
    groups_to_use: list[int] | None = None,
    background_subtract: BoolOrBoth = 'both',
    background_path: str | None = None,
    cube_build_weighting: Literal['drizzle', 'emsm', 'msm'] | None = None,
    basic_navigation: bool = False,
    step_kwargs: dict[Step, dict[str, Any]] | None = None,
    parallel_kwargs: dict[str, Any] | None = None,
    reduction_parallel_kwargs: dict[str, Any] | None = None,
    defringe: BoolOrBoth = False,
    correct_psf: BoolOrBoth = False,
    defringe_1d: BoolOrBoth = 'both',
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
        cube_build_weighting: Weighting algorithm to use when building spectral cubes.
            This is a convenience argument to set the `weighting` parameter in the
            `cube_build` step in the `stage3` pipeline, and is equivalent to passing
            a value with
            `step_kwargs={'stage3': {'steps': {'cube_build': {'weighting': ...}}}}`.
            Possible values are 'drizzle' (the default), 'emsm', and 'msm'. See
            https://jwst-pipeline.readthedocs.io/en/latest/jwst/cube_build/main.html#weighting
            for more details on the different algorithms.
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
            for the default values. The provided kwargs will be merged with the default
            kwargs.
        parallel_kwargs: Dictionary of keyword arguments to customise parallel
            processing.
        reduction_parallel_kwargs: Dictionary of keyword arguments to customise parallel
            processing for the reduction steps (i.e. `stage1`-`stage3`). This will be
            merged with `parallel_kwargs` (i.e.
            `parallel_kwargs | reduction_parallel_kwargs`).

        defringe: Toggle 2D defringing of the data. If True, defringe data will be saved
            and processed. If `'both'`, the pipeline steps will effectively be run
            twice, first without defringing, then again with defringing enabled.
            Defringed data will be saved in separate directories (e.g.
            `root_path/stage3/d1_fringe`), so both sets of data will be available for
            further analysis. Note that `defringe_1d` is generally preferred over this
            2D defringing.
        correct_psf: Toggle PSF correction of the data. If True, PSF corrected data will
            be saved and processed. If `'both'` (the default), the pipeline steps will
            effectively be run twice, first without PSF correction, then again with PSF
            correction enabled. PSF corrected data will be saved in separate directories
            (e.g. `root_path/stage3/d1_psf`), so both sets of data will be available
            for further analysis.
        defringe_1d: Toggle 1D defringing of the data. If True, 1D defringed data will
            be saved and processed. If `'both'` (the default), the pipeline steps will
            effectively be run twice, first without 1D defringing, then again with 1D
            defringing enabled. 1D defringed data will be saved in separate directories
            (e.g. `root_path/stage3/d1_fringe1d`), so both sets of data will be
            available for further analysis.
        flat_data_path: Optionally specify custom path to the flat field data. This path
            should contain `{channel}`, `{band}` and `{fringe}` placeholders, which will
            be replaced by appropriate values for each channel, band and defringe
            setting. If a directory is provided instead of a file, then the default
            flat field filename `ch{channel}-{band}_constructed_flat{fringe}.fits`
            will be appended to the path.
    """
    pipeline = MiriPipeline(
        root_path,
        parallel=parallel,
        desaturate=desaturate,
        groups_to_use=groups_to_use,
        background_subtract=background_subtract,
        background_path=background_path,
        cube_build_weighting=cube_build_weighting,
        basic_navigation=basic_navigation,
        step_kwargs=step_kwargs,
        parallel_kwargs=parallel_kwargs,
        reduction_parallel_kwargs=reduction_parallel_kwargs,
        flat_data_path=flat_data_path,
        defringe=defringe,
        defringe_1d=defringe_1d,
        correct_psf=correct_psf,
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
        defringe: BoolOrBoth = False,
        defringe_1d: BoolOrBoth = False,
        correct_psf: BoolOrBoth = 'both',
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.flat_data_path = self._make_flat_data_path(flat_data_path)
        self.defringe = defringe
        self.defringe_1d = defringe_1d
        self.correct_psf = correct_psf

    def _make_flat_data_path(self, path: str) -> str:
        """Standardise the flat data path by maybe adding a filename suffix."""
        path = self.standardise_path(path)
        if os.path.isdir(path):
            path = os.path.join(path, DEFAULT_FLAT_FILENAME)
        return path
    
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
        directories['flat'] = (
            'stage3_desaturated' if self.desaturate else 'stage3',
            directories['flat'][1],
        )
        return directories

    @property
    def stage3_file_match_hdr_keys(self) -> tuple[str, ...]:
        return ('CHANNEL', 'BAND')

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        return STAGE_DIRECTORIES_TO_PLOT

    def process_skip_steps(
        self,
        skip_steps: Collection[Step] | None,
        start_step: Step | None,
        end_step: Step | None,
    ) -> set[Step]:
        skip_steps = super().process_skip_steps(skip_steps, start_step, end_step)
        if not self.defringe:
            skip_steps.add('defringe')
        if not self.correct_psf:
            skip_steps.add('psf')
        return skip_steps

    @property
    def data_variants_individual(self) -> set[str]:
        variants = super().data_variants_individual
        if self.defringe:
            variants.add('fringe')
        if self.defringe_1d:
            variants.add('fringe1d')
        if self.correct_psf:
            variants.add('psf')
        return variants

    @property
    def data_variants_individual_required(self) -> set[str]:
        variants_required = super().data_variants_individual_required
        if self.defringe is True:
            variants_required.add('fringe')
        if self.defringe_1d is True:
            variants_required.add('fringe1d')
        if self.correct_psf is True:
            variants_required.add('psf')
        return variants_required

    def print_reduction_info(self, skip_steps: set[Step]) -> None:
        super().print_reduction_info(skip_steps)
        self.log(f'2D defringe: {self.defringe!r}', time=False)
        self.log(f'1D defringe: {self.defringe_1d!r}', time=False)
        self.log(f'Correct PSF: {self.correct_psf!r}', time=False)
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

    # defringe
    def run_defringe(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['defringe']
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(root_path, dir_in, '*cal.fits')
            if 'bg' in self.data_variants_individual:
                paths_in += self.get_paths(root_path, dir_in, 'bg', '*cal.fits')
            args_list = [(p, os.path.dirname(p), kwargs) for p in paths_in]
            runmany(
                self.defringe_fn,
                args_list,
                desc='defringe',
                **self.reduction_parallel_kwargs,
            )

    def defringe_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        ResidualFringeStep.call(
            path_in, output_dir=output_dir, save_results=True, skip=False, **kwargs
        )

    # stage3
    def get_stage3_variant_paths_in(
        self, root_path: RootPath, variant: frozenset[str]
    ) -> tuple[frozenset[str], list[str]]:
        dir_in, dir_out = self.step_directories['stage3']
        variant = variant - {'psf', 'fringe1d'}  # steps done after stage3
        if variant == frozenset():
            paths = self.get_paths(root_path, dir_in, '*cal.fits')
        elif variant == frozenset({'bg'}):
            paths = self.get_paths(root_path, dir_in, 'bg', '*cal.fits')
        elif variant == frozenset({'bg', 'fringe'}):
            paths = self.get_paths(root_path, dir_in, 'bg', '*residual_fringe.fits')
        elif variant == frozenset({'fringe'}):
            paths = self.get_paths(root_path, dir_in, '*residual_fringe.fits')
        else:
            raise ValueError(f'Unknown variant: {variant}')
        return variant, paths

    # defringe_1d
    def run_defringe_1d(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['defringe_1d']
        input_variant_combinations = {
            v - {'psf', 'fringe1d'}
            for v in self.data_variant_combinations
            if 'fringe1d' in v
        }
        for root_path in self.iterate_group_root_paths():
            paths_list: list[tuple[str, str]] = []
            for input_variant in input_variant_combinations:
                output_variant = input_variant | {'fringe1d'}
                paths_in = self.get_paths(
                    root_path,
                    dir_in,
                    '*',
                    '*_s3d.fits',
                    filter_variants=True,
                    variant_combinations={input_variant},
                )
                for p_in in paths_in:
                    # We need to add the fringe1d variant to the path
                    dirname_in = pathlib.Path(p_in).parts[-2]
                    dirname_without_variant = dirname_in.split('_')[0]
                    dirname_variants = frozenset(dirname_in.split('_')[1:])
                    dirname_out = (
                        dirname_without_variant + '_' + '_'.join(sorted(output_variant))
                    )
                    if dirname_variants != input_variant:
                        raise ValueError(
                            f'Error when checking variant path for {p_in}'
                            f'(expected variant {input_variant!r}, got {dirname_variants!r})'
                        )

                    p_out = self.replace_path_part(p_in, -3, dir_out, old=dir_in)
                    p_out = self.replace_path_part(
                        p_out, -2, dirname_out, old=dirname_in
                    )
                    if p_out == p_in:
                        raise ValueError(
                            f'Error when checking output path for {p_in}'
                            '(input and output paths are the same)'
                        )

                    paths_list.append((p_in, p_out))
            args_list = [(p_in, p_out, kwargs) for p_in, p_out in paths_list]
            runmany(
                self.defringe_1d_fn,
                args_list,
                desc='defringe_1d',
                **self.parallel_kwargs,
            )

    def defringe_1d_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        p_in, p_out, kwargs = args
        defringe_file(p_in, p_out, **kwargs)

    # stage3+spectral_leak
    def run_stage3(self, kwargs: dict[str, Any]) -> None:
        super().run_stage3(kwargs)

        # The spectral leak step needs 3A and 1B to be in the ASN, so manually run it
        # with all wavelengths in the same ASN. Putting everything in the same ASN for
        # the whole of stage3 can mess up the WCS info for moving targets, so we only
        # do it for the spectral leak step.
        try:
            spectral_leak_kwargs = kwargs['steps']['spectral_leak']
        except KeyError:
            spectral_leak_kwargs = {}
        self.run_spectral_leak(spectral_leak_kwargs)

    def run_spectral_leak(self, kwargs: dict[str, Any]) -> None:
        self.log('Running spectral leak step on any extracted 1D spectra')
        _, dir_in = self.step_directories['stage3']
        for root_path in self.iterate_group_root_paths():
            paths_list: list[tuple[str, str]] = []
            for variant in self.data_variant_combinations:
                paths_in = self.get_paths(
                    root_path,
                    dir_in,
                    '*',
                    '*_x1d.fits',
                    filter_variants=True,
                    variant_combinations={variant},
                    do_warning=False,
                )
                grouped_files = self.group_x1d_files_for_spectral_leak(paths_in)
                for (output_dir, prodname), paths_in in grouped_files.items():
                    asn_path = self.write_asn_for_stage3(
                        paths_in,
                        os.path.join(output_dir, f'{prodname}_x1d_asn.json'),
                        prodname=prodname,
                    )
                    paths_list.append((asn_path, output_dir))
            args_list = [
                (p, output_dir, kwargs) for p, output_dir in sorted(paths_list)
            ]
            if len(args_list) > 0:
                runmany(
                    self.spectral_leak_fn,
                    args_list,
                    desc='spectral_leak',
                    **self.reduction_parallel_kwargs,
                )

    def spectral_leak_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        asn_path, output_dir, kwargs = args
        container = SpectralLeakStep.call(asn_path, save_results=False, **kwargs)
        for model in container:
            # Only save the model if the spectral leak step was actually run
            if model.meta.cal_step.spectral_leak is not None:
                path = os.path.join(output_dir, model.meta.filename)
                with fits.open(path) as hdul:
                    hdr: fits.Header = hdul['PRIMARY'].header  # type: ignore
                if hdr.get('S_SPLEAK', None) == 'COMPLETE':
                    # Ensure that we don't double-apply the spectral leak step
                    print(
                        f'WARNING: {path} already has the spectral leak step applied, '
                        'skipping manual spectral leak step'
                    )
                    continue
                model.save(path)

    def group_x1d_files_for_spectral_leak(
        self, paths_in: list[str]
    ) -> dict[tuple[str, str], list[str]]:
        out: dict[tuple[str, str], list[str]] = {}
        for p_in in paths_in:
            p = pathlib.Path(p_in)
            prodname = '_'.join(p.name.split('_')[:-2])
            output_dir = str(p.parent)
            out.setdefault((output_dir, prodname), []).append(p_in)
        return out

    # psf
    def run_psf(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['psf']
        input_variant_combinations = {
            v - {'psf'} for v in self.data_variant_combinations if 'psf' in v
        }
        for root_path in self.iterate_group_root_paths():
            paths_list: list[tuple[str, str]] = []
            for input_variant in input_variant_combinations:
                output_variant = input_variant | {'psf'}
                paths_in = self.get_paths(
                    root_path,
                    dir_in,
                    '*',
                    '*_nav.fits',
                    filter_variants=True,
                    variant_combinations={input_variant},
                )
                for p_in in paths_in:
                    # We need to add the psf variant to the path
                    dirname_in = pathlib.Path(p_in).parts[-2]
                    dirname_without_variant = dirname_in.split('_')[0]
                    dirname_variants = frozenset(dirname_in.split('_')[1:])
                    dirname_out = (
                        dirname_without_variant + '_' + '_'.join(sorted(output_variant))
                    )
                    if dirname_variants != input_variant:
                        raise ValueError(
                            f'Error when checking variant path for {p_in}'
                            f'(expected variant {input_variant!r}, got {dirname_variants!r})'
                        )

                    p_out = self.replace_path_part(p_in, -3, dir_out, old=dir_in)
                    p_out = self.replace_path_part(
                        p_out, -2, dirname_out, old=dirname_in
                    )
                    if p_out == p_in:
                        raise ValueError(
                            f'Error when checking output path for {p_in}'
                            '(input and output paths are the same)'
                        )

                    paths_list.append((p_in, p_out))
            args_list = [(p_in, p_out, kwargs) for p_in, p_out in paths_list]
            runmany(
                self.psf_fn,
                args_list,
                desc='psf',
                **self.parallel_kwargs,
            )

    def psf_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        p_in, p_out, kwargs = args
        try:
            psf_correction.correct_file(p_in, p_out, **kwargs)
        except psf_correction.TargetTooSmallError:
            self.log(f'Target too small for PSF correction, skipping {p_out}')
        except psf_correction.TargetNotInFovError:
            self.log(
                f'Target not in field of view for PSF correction, skipping {p_out}'
            )

    # flat
    def run_flat(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['flat']
        self.log(f'Flat data path: {self.flat_data_path!r}', time=False)
        # Use d* to skip dither combined files
        paths_in = self.get_paths(
            self.root_path, dir_in, 'd*', '*_nav.fits', filter_variants=True
        )
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

    # plot
    def get_plot_filename_prefix(self, path: str) -> str:
        with fits.open(path) as hdul:
            hdr = hdul['PRIMARY'].header  #  type: ignore
        channel = hdr['CHANNEL']
        abc = BAND_ABC_ALIASES[hdr['BAND'].casefold().strip()]
        return f'{channel}{abc}_'


def main():
    parser = get_pipeline_argument_parser(MiriPipeline, STEP_DESCRIPTIONS, CLI_EXAMPLES)
    parser.add_argument(
        '--defringe-1d',
        action=argparse_utils.BooleanOrBothAction,
        help="""Toggle 1D defringing of the data. Note that this is generally
                preferred over the 2D defringing. If unspecified, the default is
                to run with both with and without 1D defringing.""",
    )
    parser.add_argument(
        '--defringe',
        action=argparse_utils.BooleanOrBothAction,
        help="""Toggle 2D defringing of the data. Note that 1D defringing is generally
            preferred over this 2D defringing. This 2D defringing is disabled by 
            default.""",
    )
    parser.add_argument(
        '--flat-data-path',
        type=str,
        help="""Optionally specify custom path to the flat field data. This path
            should contain `{channel}`, `{band}` and `{fringe}` placeholders, which will
            be replaced by appropriate values for each channel, band and defringe
            setting. If a directory is provided instead of a file, then the default
            flat field filename `ch{channel}-{band}_constructed_flat{fringe}.fits`
            will be appended to the path.""",
    )
    parser.add_argument(
        '--correct-psf',
        action=argparse_utils.BooleanOrBothAction,
        help="""Toggle PSF correction of the data. This PSF correction is disabled by
            default.""",
    )
    run_pipeline(**vars(parser.parse_args()))


if __name__ == '__main__':
    main()
