import argparse
import datetime
import glob
import itertools
import json
import os
import pathlib
from collections.abc import Collection
from typing import Any, Generator, Literal, TypeAlias, overload

import tqdm
from astropy.io import fits
from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline, Spec3Pipeline
from jwst.residual_fringe import ResidualFringeStep

import background_subtraction
import desaturate_data
import despike_data
import jwst_summary_animation
import jwst_summary_plots
import navigate_jwst_observations
import remove_groups
from parallel_tools import runmany
from tools import check_path

Step: TypeAlias = Literal[
    'remove_groups',
    'stage1',
    'stage2',
    'defringe',
    'stage3',
    'navigate',
    'desaturate',
    'flat',
    'despike',
    'background',
    'plot',
    'animate',
]

MIRI_STEPS = (
    'remove_groups',
    'stage1',
    'stage2',
    'defringe',
    'stage3',
    'navigate',
    'desaturate',
    'flat',
    'despike',
    'background',
    'plot',
    'animate',
)
MIRI_DEFAULT_KWARGS: dict[Step, dict[str, Any]] = {
    'stage2': {'steps': {'cube_build': {'skip': True}, 'extract_1d': {'skip': True}}},
    'defringe': {'skip': False},
    'stage3': {
        'steps': {
            'master_background': {'skip': True},
            'extract_1d': {'skip': True},
            'cube_build': {'output_type': 'band', 'coord_system': 'ifualign'},
        }
    },
}
MIRI_STAGE_DIRECTORIES_TO_PLOT = (
    'stage3',
    'stage3_desaturated',
    'stage4_flat',
    'stage5_despike',
    'stage6_background',
)

NIRSPEC_STEPS = (
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
)
NIRSPEC_DEFAULT_KWARGS: dict[Step, dict[str, Any]] = {
    'stage2': {
        'steps': {
            'cube_build': {'coord_system': 'ifualign'},
            'extract_1d': {'skip': True},
        }
    },
    'stage3': {
        'steps': {
            'cube_build': {'coord_system': 'ifualign'},
            'extract_1d': {'skip': True},
        }
    },
}
NIRSPEC_STAGE_DIRECTORIES_TO_PLOT = (
    'stage3',
    'stage3_desaturated',
    'stage4_despike',
    'stage5_background',
)


class Pipeline:
    def __init__(
        self,
        root_path: str,
        *,
        parallel: float | bool = False,
        desaturate: bool = True,
        groups_to_use: list[int] | None = None,
        background_path: str | None = None,
        basic_navigation: bool = False,
        step_kwargs: dict[Step, dict[str, Any]] | None = None,
        parallel_kwargs: dict[str, Any] | None = None,
        reduction_parallel_kwargs: dict[str, Any] | None = None,
    ):
        # TODO docstring
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
        self.parallel_kwargs = parallel_kwargs
        self.reduction_parallel_kwargs = reduction_parallel_kwargs
        self.step_kwargs = step_kwargs or {}
        self.root_path = self.standardise_path(root_path)
        self.desaturate = desaturate
        self.groups_to_use = groups_to_use
        self.background_path = self.standardise_path(background_path)
        self.basic_navigation = basic_navigation

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.root_path!r})'

    # Pipeline data to override in subclasses ------------------------------------------
    @property
    def instrument(self) -> str:
        raise NotImplementedError

    @property
    def steps(self) -> tuple[Step, ...]:
        raise NotImplementedError

    @property
    def default_kwargs(self) -> dict[Step, dict[str, Any]]:
        raise NotImplementedError

    @property
    def step_directories(self) -> dict[Step, tuple[str, str]]:
        """
        Dictionary containing input and output directories for each step.
        """
        return {
            'remove_groups': ('stage0', 'stage0'),
            'stage1': ('stage0', 'stage1'),
            'stage2': ('stage1', 'stage2'),
            'defringe': ('stage2', 'stage2'),
            'stage3': ('stage2', 'stage3'),
            'navigate': ('stage3', 'stage3'),
            'desaturate': ('stage3', 'stage3_desaturated'),
            # later steps overriden in subclasses
        }

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        raise NotImplementedError

    # Pipeline running methods ---------------------------------------------------------
    def run(
        self,
        *,
        skip_steps: Collection[Step] | None = None,
        start_step: Step | None = None,
        end_step: Step | None = None,
    ):
        # TODO docstring
        skip_steps = self.process_skip_steps(skip_steps, start_step, end_step)
        self.log(f'Running {self.instrument} pipeline')
        self.print_reduction_info(skip_steps)
        print()
        for step in self.steps:
            if step in skip_steps:
                continue
            self.run_step(step)

        self.log(f'{self.instrument} pipeline completed with')
        self.print_reduction_info(skip_steps)
        self.log('ALL DONE')

    def run_step(self, step: Step) -> None:
        """Run pipeline step"""
        self.log(f'Running {step} step')
        kwargs = self.default_kwargs.get(step, {}) | self.step_kwargs.get(step, {})
        if kwargs:
            self.log(f'Arguments: {kwargs!r}', time=False)
        getattr(self, f'run_{step}')(kwargs)
        self.log(f'{step} step completed')

    def process_skip_steps(
        self,
        skip_steps: Collection[Step] | None,
        start_step: Step | None,
        end_step: Step | None,
    ) -> set[Step]:
        """
        Process the various step arguments to create a comprehensive set of steps to
        skip.
        """
        if skip_steps is None:
            skip_steps = []
        for s in skip_steps:
            if s not in self.steps:
                raise ValueError(
                    f'Invalid skip step: {s!r} (valid steps: {self.steps})'
                )
        skip_steps = set(skip_steps)
        if start_step is not None:
            if start_step not in self.steps:
                raise ValueError(
                    f'Invalid start step: {start_step!r} (valid steps: {self.steps})'
                )
            skip_steps.update(self.steps[: self.steps.index(start_step)])
        if end_step is not None:
            if end_step not in self.steps:
                raise ValueError(
                    f'Invalid end step: {end_step!r} (valid steps: {self.steps})'
                )
            skip_steps.update(self.steps[self.steps.index(end_step) + 1 :])

        # Add skip steps dependent on provided arguments
        if not self.desaturate:
            skip_steps.add('desaturate')
        if self.background_path is None:
            skip_steps.add('background')  # TODO update this if the bg step changes?
        if self.basic_navigation:
            skip_steps.add('animate')
        return skip_steps

    def print_reduction_info(self, skip_steps: set[Step]) -> None:
        """Print information about the current reduction run."""
        self.log(f'Root path: {self.root_path!r}', time=False)
        if skip_steps:
            self.log(
                f'Skipping steps: {", ".join([repr(s) for s in self.steps if s in skip_steps])}',
                time=False,
            )
        else:
            self.log('Running all pipeline steps', time=False)
        self.log(f'Desaturate: {self.desaturate!r}', time=False)
        if self.groups_to_use:
            self.log(f'Groups to keep: {self.groups_to_use!r}', time=False)
        if self.basic_navigation:
            self.log(f'Basic navigation: {self.basic_navigation!r}', time=False)

    # Utility methods ------------------------------------------------------------------
    @staticmethod
    def log(*messages: Any, time: bool = True) -> None:
        """
        Print a message with a timestamp.

        Args:
            *messages: Messages passed to `print()`.
            time: Toggle showing the timestamp.
        """
        prefix = datetime.datetime.now().strftime('%H:%M:%S') if time else ' ' * 8
        print(prefix, *messages, flush=True)

    # path processing...
    @staticmethod
    @overload
    def standardise_path(path: str) -> str:
        ...

    @staticmethod
    @overload
    def standardise_path(path: None) -> None:
        ...

    @staticmethod
    def standardise_path(path: str | None) -> str | None:
        """Standardise a path by expanding environment variables and user."""
        if path is None:
            return None
        return os.path.expandvars(os.path.expanduser(path))

    @staticmethod
    def replace_path_part(
        path: str, idx: int, new: str, *, old: str | None = None
    ) -> str:
        """
        Replace a part of a path with a new value.

        Args:
            path: The path to modify.
            idx: The index of the part to replace in the directory tree.
            new: The new value to use.
            old: If not None, check that the value at the index is this value before
                replacing it.

        Returns:
            The modified path.
        """
        p = pathlib.Path(path)
        parts = list(p.parts)
        if old is not None and parts[idx] != old:
            raise ValueError(f'Expected {old!r} at index {idx} in {path!r}')
        parts[idx] = new
        return str(pathlib.Path(*parts))

    @staticmethod
    def replace_path_suffix(path: str, new: str, *, old: str | None = None) -> str:
        """
        Replace the suffix of a path with a new value.

        Args:
            path: The path to modify.
            new: The new value to use.
            old: If not None, check that the suffix is this value before replacing it.

        Returns:
            The modified path.
        """
        p = pathlib.Path(path)
        if old is not None and p.suffix != old:
            raise ValueError(f'Expected {old!r} suffix in {path!r}')
        return str(p.with_suffix(new))

    # path getting/filtering...
    def get_paths(self, *path_parts: str, filter_variants: bool = False) -> list[str]:
        """Get a list of paths matching the given path parts."""
        paths = sorted(glob.glob(os.path.join(self.root_path, *path_parts)))
        if filter_variants:
            paths = self.filter_paths_for_data_variants(paths)
        return paths

    @property
    def group_root_paths(self) -> list[str]:
        group_root_paths = [self.root_path]
        if self.desaturate:
            # If desaturating, we also need to reduce the data with fewer groups
            # This list is sorted such that the number of groups is decreasing (so that the
            # desaturation works correctly)
            reduced_group_root_paths = sorted(
                glob.glob(os.path.join(self.root_path, 'groups', '*_groups')),
                reverse=True,
                key=lambda _p: int(os.path.basename(_p).split('_')[0]),
            )
            if self.groups_to_use is not None:
                reduced_group_root_paths = [
                    _p
                    for _p in reduced_group_root_paths
                    if int(os.path.basename(_p).split('_')[0]) in self.groups_to_use
                ]
            group_root_paths.extend(reduced_group_root_paths)
        return group_root_paths

    def iterate_group_root_paths(self) -> Generator[str, Any, None]:
        group_root_paths = self.group_root_paths
        for idx, root_path in enumerate(group_root_paths):
            if len(group_root_paths) > 1:
                self.log(
                    f'Processing group directory {idx+1}/{len(group_root_paths)} ({root_path!r})'
                )
            yield root_path

    @property
    def data_variants_individual(self) -> set[str]:
        """
        Possible data variants that are applied to the dithered folders.

        E.g. if the individual data variants are {'bg', 'fringe'}, then the following
        the dithered folders will be created:
            - d1
            - d1_bg
            - d1_fringe
            - d1_bg_fringe
        etc.
        """
        return set()

    @property
    def data_variant_combinations(self) -> set[frozenset[str]]:
        """
        All possible combinations of data variants.
        """
        # Start off with empty set (i.e. no variants)
        combinations: set[frozenset[str]] = {frozenset()}
        for variant in self.data_variants_individual:
            # For each variant, add it to each existing combination
            combinations.update({c | {variant} for c in combinations})
        return combinations

    def test_path_for_data_variants(self, path: str) -> bool:
        """
        Test if the given path contains the valid data variants.

        Args:
            path: Path to test.

        Returns:
            True if the path contains the valid data variants, False otherwise.
        """
        dirname = pathlib.Path(path).parts[-2]
        # dirnames are of format d1_bg_fringe, so split on _ to get variants and drop
        # the first element (the dither number/combination). Use a set to ensure order
        # doesn't matter
        path_variants = frozenset(dirname.split('_')[1:])
        return path_variants in self.data_variant_combinations

    def filter_paths_for_data_variants(self, paths: list[str]) -> list[str]:
        """
        Filter a list of paths to only include those containing the valid data variants.

        Args:
            paths: List of paths to filter.

        Returns:
            Filtered list of paths.
        """
        return [p for p in paths if self.test_path_for_data_variants(p)]

    # Pipeline steps -------------------------------------------------------------------
    # remove_groups
    def run_remove_groups(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['remove_groups']
        paths_in = self.get_paths(dir_in, '*_uncal.fits')
        self.log(f'Processing {len(paths_in)} files...', time=False)
        for p in tqdm.tqdm(paths_in, desc='remove_groups'):
            remove_groups.remove_groups_from_file(p, self.groups_to_use)

    # stage1
    def run_stage1(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['remove_groups']
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(dir_in, '*.fits')
            output_dir = os.path.join(root_path, dir_out)
            args_list = [(p, output_dir, kwargs) for p in paths_in]
            self.log(f'Output directory: {output_dir!r}', time=False)
            check_path(output_dir)
            runmany(
                self.reduction_detector1_fn,
                args_list,
                desc='stage1',
                **self.reduction_parallel_kwargs,
            )

    @staticmethod
    def reduction_detector1_fn(args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        Detector1Pipeline.call(
            path_in, output_dir=output_dir, save_results=True, **kwargs
        )

    # stage2
    def run_stage2(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['stage2']
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(dir_in, '*rate.fits')
            output_dir = os.path.join(root_path, dir_out)
            args_list = [(p, output_dir, kwargs) for p in paths_in]
            self.log(f'Output directory: {output_dir!r}', time=False)
            check_path(output_dir)
            runmany(
                self.reduction_spec2_fn,
                args_list,
                desc='stage2',
                **self.reduction_parallel_kwargs,
            )

    @staticmethod
    def reduction_spec2_fn(args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        Spec2Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)

    # TODO background

    # stage3
    def run_stage3(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['stage3']
        for root_path in self.iterate_group_root_paths():
            # TODO deal with variants here
            # TODO stage3 stuff
            pass

    # navigate
    def run_navigate(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['navigate']
        self.log(f'Basic navigation: {self.basic_navigation}')
        navigate_jwst_observations.load_kernels()
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(
                root_path, dir_in, '*', '*_s3d.fits', filter_variants=True
            )
            navigate_jwst_observations.navigate_multiple(
                *paths_in, basic=self.basic_navigation, rename_directory=False, **kwargs
            )

    # desaturate
    def run_desaturate(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['desaturate']

        # Group paths by their relative path to the group root path
        paths_dict: dict[str, list[str]] = {}
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(
                root_path, dir_in, '*', '*_nav.fits', filter_variants=True
            )
            for p_in in paths_in:
                # get file path relative to the root of the group root path
                relpath = os.path.relpath(p_in, root_path)
                relpath = self.replace_path_part(relpath, -3, dir_out, old=dir_in)
                p_out = os.path.join(root_path, relpath)
                paths_dict.setdefault(p_out, []).append(p_in)

        args_list = [
            (paths_in, p_out, kwargs) for p_out, paths_in in paths_dict.items()
        ]
        runmany(
            self.desaturate_fn, args_list, desc='desaturate', **self.parallel_kwargs
        )

    @staticmethod
    def desaturate_fn(args: tuple[list[str], str, dict[str, Any]]) -> None:
        paths_in, path_out, kwargs = args
        desaturate_data.replace_saturated(paths_in, path_out, **kwargs)

    # despike
    def run_despike(self, kwargs: dict[str, Any]) -> None:
        if self.parallel_kwargs.get('parallel_frac', False):
            # don't want progress bars to be printed when running in parallel
            kwargs.setdefault('progress_bar', False)

        dir_in, dir_out = self.step_directories['despike']
        paths_in = self.get_paths(
            self.root_path, dir_in, '*', '*_nav.fits', filter_variants=True
        )
        paths_out = [
            self.replace_path_part(p, -3, dir_out, old=dir_in) for p in paths_in
        ]
        args_list = [(p_in, p_out, kwargs) for p_in, p_out in zip(paths_in, paths_out)]
        runmany(self.despike_fn, args_list, desc='despike', **self.parallel_kwargs)

    @staticmethod
    def despike_fn(args: tuple[str, str, dict[str, Any]]) -> None:
        p_in, p_out, kwargs = args
        despike_data.despike_cube(p_in, p_out, **kwargs)

    # plot
    def run_plot(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['plot']
        for root_path in self.iterate_group_root_paths():
            # TODO deal with stages to plot here
            # TODO add prefixes here
            pass
        
    # TODO animate

class MiriPipeline(Pipeline):
    @property
    def instrument(self) -> str:
        return 'MIRI'

    @property
    def steps(self) -> tuple[Step, ...]:
        return MIRI_STEPS

    @property
    def default_kwargs(self) -> dict[Step, dict[str, Any]]:
        return MIRI_DEFAULT_KWARGS

    @property
    def step_directories(self) -> dict[Step, tuple[str, str]]:
        return super().step_directories | {}  # TODO

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        return MIRI_STAGE_DIRECTORIES_TO_PLOT

    def process_skip_steps(
        self,
        skip_steps: Collection[Step] | None,
        start_step: Step | None,
        end_step: Step | None,
    ) -> set[Step]:
        skip_steps = super().process_skip_steps(skip_steps, start_step, end_step)
        # TODO defringe stuff
        return skip_steps

    @property
    def data_variants_individual(self) -> set[str]:
        return super().data_variants_individual | {'fringe'}

    # Step overrides
    @staticmethod
    def reduction_detector1_fn(args: tuple[str, str, dict[str, Any]]) -> None:
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
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths('stage2', '*cal.fits')
            output_dir = os.path.join(root_path, 'stage2')
            args_list = [(p, output_dir, kwargs) for p in paths_in]
            self.log(f'Output directory: {output_dir!r}', time=False)
            self.log(f'Processing {len(args_list)} files...', time=False)
            check_path(output_dir)
            runmany(
                self.defringe_fn,
                args_list,
                desc='defringe',
                **self.reduction_parallel_kwargs,
            )

    @staticmethod
    def defringe_fn(args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        ResidualFringeStep.call(
            path_in, output_dir=output_dir, save_results=True, **kwargs
        )


class NirspecPipeline(Pipeline):
    @property
    def instrument(self) -> str:
        return 'NIRSpec'

    @property
    def steps(self) -> tuple[Step, ...]:
        return NIRSPEC_STEPS

    @property
    def default_kwargs(self) -> dict[Step, dict[str, Any]]:
        return NIRSPEC_DEFAULT_KWARGS

    @property
    def step_directories(self) -> dict[Step, tuple[str, str]]:
        return super().step_directories | {}  # TODO

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        return NIRSPEC_STAGE_DIRECTORIES_TO_PLOT

    # Step overrides
    @staticmethod
    def reduction_spec2_fn(args: tuple[str, str, dict[str, Any]]) -> None:
        try:
            return super().reduction_spec2_fn(args)
        except NoDataOnDetectorError:
            path_in, output_dir, kwargs = args
            print(f'No data on detector for {path_in!r}, skipping')
