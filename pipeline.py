import argparse
import datetime
import glob
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


class Pipeline:
    STEPS: tuple[Step, ...] = ()  # override this in subclasses
    DEFAULT_KWARGS: dict[Step, dict[str, Any]] = {}  # override this in subclasses
    DATA_STAGES: tuple[str, ...] = ()  # override this in subclasses

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

    @property
    def instrument(self) -> str:
        raise NotImplementedError

    # Pipeline running
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
        for step in self.STEPS:
            if step in skip_steps:
                continue
            self.run_step(step)

        self.log(f'{self.instrument} pipeline completed with')
        self.print_reduction_info(skip_steps)
        self.log('ALL DONE')

    def run_step(self, step: Step) -> None:
        """Run pipeline step"""
        self.log(f'Running {step} step')
        kwargs = self.DEFAULT_KWARGS.get(step, {}) | self.step_kwargs.get(step, {})
        if kwargs:
            self.log(f'Arguments: {kwargs!r}', time=False)
        skipped = getattr(self, f'run_{step}')(kwargs)
        if not skipped:
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
            if s not in self.STEPS:
                raise ValueError(
                    f'Invalid skip step: {s!r} (valid steps: {self.STEPS})'
                )
        skip_steps = set(skip_steps)
        if start_step is not None:
            if start_step not in self.STEPS:
                raise ValueError(
                    f'Invalid start step: {start_step!r} (valid steps: {self.STEPS})'
                )
            skip_steps.update(self.STEPS[: self.STEPS.index(start_step)])
        if end_step is not None:
            if end_step not in self.STEPS:
                raise ValueError(
                    f'Invalid end step: {end_step!r} (valid steps: {self.STEPS})'
                )
            skip_steps.update(self.STEPS[self.STEPS.index(end_step) + 1 :])
        if self.basic_navigation:
            skip_steps.add('animate')
        return skip_steps

    def print_reduction_info(self, skip_steps: set[Step]) -> None:
        """Print information about the current reduction run."""
        self.log(f'Root path: {self.root_path!r}', time=False)
        if skip_steps:
            self.log(
                f'Skipping steps: {", ".join([repr(s) for s in self.STEPS if s in skip_steps])}',
                time=False,
            )
        else:
            self.log('Running all pipeline steps', time=False)
        self.log(f'Desaturate: {self.desaturate!r}', time=False)
        if self.groups_to_use:
            self.log(f'Groups to keep: {self.groups_to_use!r}', time=False)
        if self.basic_navigation:
            self.log(f'Basic navigation: {self.basic_navigation!r}', time=False)

    # Utility methods
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
    def log(*messages: Any, time: bool = True) -> None:
        """
        Print a message with a timestamp.

        Args:
            *messages: Messages passed to `print()`.
            time: Toggle showing the timestamp.
        """
        prefix = datetime.datetime.now().strftime('%H:%M:%S') if time else ' ' * 8
        print(prefix, *messages, flush=True)

    def get_paths(self, *path_parts: str) -> list[str]:
        """Get a list of paths matching the given path parts."""
        return sorted(glob.glob(os.path.join(self.root_path, *path_parts)))

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
    def data_variants(self) -> list[str]:
        # TODO some sort of list of variants of the data e.g.
        # '', _fringe, _bg, _bg_fringe
        pass

    def iterate_data_variants(self) -> Generator[str, Any, None]:
        pass

    # Pipeline steps
    def run_remove_groups(self, kwargs: dict[str, Any]) -> None:
        paths_in = self.get_paths('stage0', '*_uncal.fits')
        self.log(f'Processing {len(paths_in)} files...', time=False)
        for p in tqdm.tqdm(paths_in, desc='remove_groups'):
            remove_groups.remove_groups_from_file(p, self.groups_to_use)

    def run_stage1(self, kwargs: dict[str, Any]) -> None:
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths('stage0', '*.fits')
            output_dir = os.path.join(root_path, 'stage1')
            args_list = [(p, output_dir, kwargs) for p in paths_in]
            self.log(f'Output directory: {output_dir!r}', time=False)
            self.log(f'Processing {len(args_list)} files...', time=False)
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

    def run_stage2(self, kwargs: dict[str, Any]) -> None:
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths('stage1', '*rate.fits')
            output_dir = os.path.join(root_path, 'stage2')
            args_list = [(p, output_dir, kwargs) for p in paths_in]
            self.log(f'Output directory: {output_dir!r}', time=False)
            self.log(f'Processing {len(args_list)} files...', time=False)
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


class MIRIPipeline(Pipeline):
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
        'background',
        'plot',
        'animate',
    )
    DEFAULT_KWARGS = {
        'stage2': {
            'steps': {'cube_build': {'skip': True}, 'extract_1d': {'skip': True}}
        },
        'defringe': {'skip': False},
        'stage3': {
            'steps': {
                'master_background': {'skip': True},
                'extract_1d': {'skip': True},
                'cube_build': {'output_type': 'band', 'coord_system': 'ifualign'},
            }
        },
    }
    DATA_STAGES = (
        'stage3',
        'stage3_desaturated',
        'stage4_flat',
        'stage5_despike',
        'stage6_background',
    )

    @property
    def instrument(self) -> str:
        return 'MIRI'

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


class NIRSpecPipeline(Pipeline):
    STEPS = (
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
    DEFAULT_KWARGS = {
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
    DATA_STAGES = (
        'stage3',
        'stage3_desaturated',
        'stage4_despike',
        'stage5_background',
    )

    @property
    def instrument(self) -> str:
        return 'NIRSpec'

    @staticmethod
    def reduction_spec2_fn(args: tuple[str, str, dict[str, Any]]) -> None:
        try:
            return super().reduction_spec2_fn(args)
        except NoDataOnDetectorError:
            path_in, output_dir, kwargs = args
            print(f'No data on detector for {path_in!r}, skipping')
