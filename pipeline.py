"""
Common JWST pipeline code used by the MIRI and NIRSpec pipelines.

GitHub repository: https://github.com/JWSTGiantPlanets/pipelines
"""
import argparse
import datetime
import glob
import json
import os
import pathlib
from typing import (
    Any,
    Collection,
    Generator,
    Literal,
    NewType,
    Type,
    TypeAlias,
    cast,
    overload,
)

import tqdm
from astropy.io import fits
from jwst.associations.asn_from_list import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline, Spec3Pipeline

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
    'stage3',
    'navigate',
    'desaturate',
    'despike',
    'plot',
    'animate',
    'psf',  # MIRI only
    'flat',  # MIRI only
    'defringe',  # MIRI only
]
BoolOrBoth: TypeAlias = Literal[True, False, 'both']
RootPath = NewType('RootPath', str)


STEP_DIRECTORIES: dict[Step, tuple[str, str]] = {
    'remove_groups': ('stage0', 'stage0'),
    'stage1': ('stage0', 'stage1'),
    'stage2': ('stage1', 'stage2'),
    'stage3': ('stage2', 'stage3'),
    'navigate': ('stage3', 'stage3'),
    'desaturate': ('stage3', 'stage3_desaturated'),
    # later steps overriden in subclasses
    'plot': ('', 'plots'),
    'animate': ('', 'animation'),
}


class Pipeline:
    """
    Base class for JWST pipeline objects.

    This cannot be run directly, but should be subclassed for each instrument.

    Args:
        root_path: Path to the root directory containing the data. This directory should
            contain a subdirectory with the stage 0 data (i.e. `root_path/stage0`).
            Additional directories will be created automatically for each step of the
            pipeline (e.g. `root_path/stage1`, `root_path/plots` etc.).
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
            `stage2`.
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
            can will customise the `stage3` step.
        parallel_kwargs: Dictionary of keyword arguments to customise parallel
            processing.
        reduction_parallel_kwargs: Dictionary of keyword arguments to customise parallel
            processing for the reduction steps (i.e. `stage1`-`stage3`). This will be
            merged with `parallel_kwargs` (i.e.
            `parallel_kwargs | reduction_parallel_kwargs`).
    """

    def __init__(
        self,
        root_path: str,
        *,
        parallel: float | bool = False,
        desaturate: bool = True,
        groups_to_use: list[int] | None | str = None,
        background_subtract: BoolOrBoth,
        background_path: str | None = None,
        basic_navigation: bool = False,
        step_kwargs: dict[Step, dict[str, Any]] | None | str = None,
        parallel_kwargs: dict[str, Any] | None = None,
        reduction_parallel_kwargs: dict[str, Any] | None = None,
    ):
        # Process CLI arguements
        if isinstance(groups_to_use, str):
            groups_to_use = [int(g) for g in groups_to_use.split(',')]
        if isinstance(step_kwargs, str):
            step_kwargs = cast(dict[Step, dict[str, Any]], json.loads(step_kwargs))

        parallel_kwargs = dict(
            parallel_frac=parallel,
            timeout=2 * 60 * 60,  # avgerage of 2 hours/job
        ) | (parallel_kwargs or {})
        reduction_parallel_kwargs = (
            parallel_kwargs
            | dict(
                timeout=5 * 60 * 60,  # average of 5 hours/job
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
        self.root_path = RootPath(self.standardise_path(root_path))
        self.desaturate = desaturate
        self.groups_to_use = groups_to_use
        self.background_subtract = background_subtract
        self.background_path = self.standardise_path(background_path)
        self.basic_navigation = basic_navigation

        for k in self.step_kwargs.keys():
            if k not in self.steps:
                raise ValueError(f'step_kwargs contains invalid step name {k!r}')

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.root_path!r})'

    # Pipeline data to override in subclasses ------------------------------------------
    @staticmethod
    def get_instrument() -> str:
        """String giving the instrument name."""
        raise NotImplementedError

    @property
    def steps(self) -> tuple[Step, ...]:
        """Tuple of step names in the order that they should be executed."""
        raise NotImplementedError

    @property
    def default_kwargs(self) -> dict[Step, dict[str, Any]]:
        """
        Dictionary of default keyword arguments for each step. These will be merged
        with the custom `step_kwargs` argument passed to the pipeline.
        """
        raise NotImplementedError

    @property
    def step_directories(self) -> dict[Step, tuple[str, str]]:
        """
        Dictionary containing input and output directories for each step.
        """
        return STEP_DIRECTORIES.copy()

    @property
    def stage3_file_match_hdr_keys(self) -> tuple[str, ...]:
        """Header keys for matching stage2 files to combine in stage3."""
        raise NotImplementedError

    @property
    def stage_directories_to_plot(self) -> tuple[str, ...]:
        """Tuple of directory names to generate plots for."""
        raise NotImplementedError

    # Pipeline running methods ---------------------------------------------------------
    def run(
        self,
        *,
        skip_steps: Collection[Step] | None = None,
        start_step: Step | None = None,
        end_step: Step | None = None,
    ):
        """
        Run the pipeline.

        Args:
            skip_steps: List of steps to skip.
            start_step: Convenience argument to add all steps before `start_step` to
                `skip_steps`.
            end_step: Convenience argument to add all steps after `end_step` to
                `skip_steps`.
        """
        skip_steps = self.process_skip_steps(skip_steps, start_step, end_step)
        self.log(f'Running {self.get_instrument()} pipeline')
        self.print_reduction_info(skip_steps)

        for step in self.steps:
            if step in skip_steps:
                continue
            self.run_step(step)

        print()
        self.log(f'{self.get_instrument()} pipeline completed with')
        self.print_reduction_info(skip_steps)
        self.log('ALL DONE')

    def run_step(self, step: Step) -> None:
        """Run individual pipeline step."""
        print()
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
        Process the various arguments to create a comprehensive set of steps to skip.
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
            skip_steps.add('remove_groups')
            skip_steps.add('desaturate')
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
        self.log(f'Background subtract: {self.background_subtract!r}', time=False)
        self.log(f'Background path: {self.background_path!r}', time=False)
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
    def get_paths(
        self,
        root: RootPath,
        *path_parts: str,
        filter_variants: bool = False,
        variant_combinations: set[frozenset[str]] | None = None,
    ) -> list[str]:
        """Get a list of paths matching the given path parts."""
        paths = sorted(glob.glob(os.path.join(root, *path_parts)))
        if filter_variants:
            paths = self.filter_paths_for_data_variants(paths, variant_combinations)
        return paths

    @property
    def group_root_paths(self) -> list[RootPath]:
        """
        List of relative root paths for different numbers of reduced groups, in
        descending order.

        If desaturation is enabled, this list will be of the form:
        `[root_path, root_path+'/groups/4_groups', root_path+'/groups/3_groups', ...]`.
        If desaturation is disabled, then the returned list will be `[root_path]`.
        """
        group_root_paths: list[str] = [self.root_path]
        if self.desaturate:
            # If desaturating, we also need to reduce the data with fewer groups
            # This list is sorted such that the number of groups is decreasing (so that
            # the desaturation works correctly)
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
        return group_root_paths  # type: ignore

    def iterate_group_root_paths(self) -> Generator[RootPath, Any, None]:
        """
        Iterate over the group root paths, yielding each path and printing a message.
        """
        group_root_paths = self.group_root_paths
        for idx, root_path in enumerate(group_root_paths):
            if len(group_root_paths) > 1:
                self.log(
                    f'Processing group directory {idx+1}/{len(group_root_paths)}: {root_path!r}'
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

        See also the `data_variants_individual_required` property.
        """
        variants = set()
        if self.background_subtract:
            variants.add('bg')
        return variants

    @property
    def data_variants_individual_required(self) -> set[str]:
        """
        Data variants that must be included in the dithered folders.

        This must be a subset of the `data_variants_individual` property.

        E.g. if the individual data variants are {'bg', 'fringe'}, and the required
        variants are {'bg'}, then the following the dithered folders will be created:
            - d1_bg
            - d1_bg_fringe
        i.e. any variant which does not contain 'bg' will not be created.
        """
        variants_required = set()
        if self.background_subtract is True:
            variants_required.add('bg')
        return variants_required

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

        # Remove any combinations that don't include the required variants
        combinations = {
            c for c in combinations if self.data_variants_individual_required <= c
        }

        return combinations

    def test_path_for_data_variants(
        self, path: str, variant_combinations: set[frozenset[str]] | None = None
    ) -> bool:
        """
        Test if the given path contains the valid data variants.

        Args:
            path: Path to test.
            variant_combinations: Set of variant combinations to test against. If None,
                then the `data_variant_combinations` property is used.

        Returns:
            True if the path contains the valid data variants, False otherwise.
        """
        if variant_combinations is None:
            variant_combinations = self.data_variant_combinations
        dirname = pathlib.Path(path).parts[-2]
        # dirnames are of format d1_bg_fringe, so split on _ to get variants and drop
        # the first element (the dither number/combination). Use a set to ensure order
        # doesn't matter
        path_variants = frozenset(dirname.split('_')[1:])
        return path_variants in variant_combinations

    def filter_paths_for_data_variants(
        self, paths: list[str], variant_combinations: set[frozenset[str]] | None = None
    ) -> list[str]:
        """
        Filter a list of paths to only include those containing the valid data variants.

        Args:
            paths: List of paths to filter.
            variant_combinations: Set of variant combinations to test against. If None,
                then the `data_variant_combinations` property is used.

        Returns:
            Filtered list of paths.
        """
        return [
            p
            for p in paths
            if self.test_path_for_data_variants(p, variant_combinations)
        ]

    def get_file_match_key(
        self, path: str | fits.Header, hdr_keys: tuple[str, ...]
    ) -> tuple[str | None, ...]:
        """
        Get a key used to e.g. match background and science files.
        """
        if isinstance(path, str):
            with fits.open(path) as hdul:
                hdr = hdul['PRIMARY'].header  #  type: ignore
        else:
            hdr = path
        return tuple(hdr.get(k, None) for k in hdr_keys)  # type: ignore

    # Pipeline steps -------------------------------------------------------------------
    # remove_groups
    # pylint: disable-next=unused-argument
    def run_remove_groups(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['remove_groups']
        paths_in = self.get_paths(self.root_path, dir_in, '*_uncal.fits')
        self.log(f'Processing {len(paths_in)} files...', time=False)
        for p in tqdm.tqdm(paths_in, desc='remove_groups'):
            remove_groups.remove_groups_from_file(p, self.groups_to_use)

    # stage1
    def run_stage1(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['stage1']
        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(root_path, dir_in, '*uncal.fits')
            output_dir = os.path.join(root_path, dir_out)
            args_list = [(p, output_dir, kwargs) for p in paths_in]
            check_path(output_dir)
            runmany(
                self.reduction_detector1_fn,
                args_list,
                desc='stage1',
                **self.reduction_parallel_kwargs,
            )

    def reduction_detector1_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        Detector1Pipeline.call(
            path_in, output_dir=output_dir, save_results=True, **kwargs
        )

    # stage2
    def run_stage2(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['stage2']

        if self.background_subtract and self.background_path is None:
            self.log(
                'No background path provided so skipping background subtraction',
                time=False,
            )
            background_options = [False]
        elif self.background_subtract == 'both':
            self.log(
                'Producing background subtracted and non-background subtracted data',
                time=False,
            )
            background_options = [True, False]
        elif self.background_subtract == True:
            self.log('Only producing background subtracted data', time=False)
            background_options = [True]
        else:
            self.log('Only producing non-background subtracted data', time=False)
            background_options = [False]

        if True in background_options:
            self.log(f'Background path: {self.background_path!r}', time=False)
            background_path_dict = self.get_background_path_dict()
            self.log(
                f'Found {len(background_path_dict)} background file groups to consider',
            )
        else:
            background_path_dict = {}

        for root_path in self.iterate_group_root_paths():
            paths_in = self.get_paths(root_path, dir_in, '*rate.fits')
            args_list: list[tuple[str, str, dict[str, Any]]] = []
            for background in background_options:
                output_dir = os.path.join(root_path, dir_out)
                if background:
                    output_dir = os.path.join(output_dir, 'bg')
                check_path(output_dir)
                for p_in in paths_in:
                    if background:
                        try:
                            bg_paths = background_path_dict[
                                self.get_file_match_key_for_background(p_in)
                            ]
                        except KeyError:
                            self.log(
                                f'WARNING: no background file found for {p_in!r}, skipping'
                            )
                            continue
                    else:
                        bg_paths = None
                    asn_path = self.write_asn_for_stage2(p_in, bg_paths)
                    args_list.append((asn_path, output_dir, kwargs))

            runmany(
                self.reduction_spec2_fn,
                args_list,
                desc='stage2',
                **self.reduction_parallel_kwargs,
            )

    def reduction_spec2_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        path_in, output_dir, kwargs = args
        Spec2Pipeline.call(path_in, output_dir=output_dir, save_results=True, **kwargs)

    def get_background_path_dict(self) -> dict[tuple[str | None, ...], list[str]]:
        """
        Get dictionary of background paths keyed by the match key.
        """
        if self.background_path is None:
            return {}
        dir_in, dir_out = self.step_directories['stage2']
        paths = sorted(
            glob.glob(os.path.join(self.background_path, dir_in, '*rate.fits'))
        )
        out: dict[tuple[str | None, ...], list[str]] = {}
        for p in paths:
            key = self.get_file_match_key_for_background(p)
            out.setdefault(key, []).append(p)
        return out

    def get_file_match_key_for_background(
        self, path: str | fits.Header
    ) -> tuple[str | None, ...]:
        """
        Get a key used to match background and science files.
        """
        return self.get_file_match_key(
            path, ('INSTRUME', 'CHANNEL', 'BAND', 'DETECTOR', 'FILTER', 'GRATING')
        )

    def write_asn_for_stage2(
        self,
        science_path: str,
        bg_paths: list[str] | None = None,
    ) -> str:
        prodname = os.path.basename(science_path).replace('_rate.fits', '')
        asn = asn_from_list([science_path], product_name=prodname)
        if bg_paths is not None:
            for p_bg in bg_paths:
                asn['products'][0]['members'].append(
                    {'expname': p_bg, 'exptype': 'background'}
                )
        asn_path = os.path.join(
            os.path.dirname(science_path),
            os.path.basename(science_path).replace(
                '.fits', ('_bg' if bg_paths is not None else '') + '_asn.json'
            ),
        )
        _, serialized = asn.dump()
        with open(asn_path, 'w', encoding='utf-8') as f:
            f.write(serialized)
        return asn_path

    # stage3
    def run_stage3(self, kwargs: dict[str, Any]) -> None:
        dir_in, dir_out = self.step_directories['stage3']
        for root_path in self.iterate_group_root_paths():
            paths_list: list[tuple[str, str]] = []
            variants_done = set()
            for full_variant in self.data_variant_combinations:
                variant, paths_in = self.get_stage3_variant_paths_in(
                    root_path, full_variant
                )
                # Skip if we have already done this variant. E.g. for MIRI, the PSF
                # variants occur after stage3, so fulll_variant={'bg'} and
                # full_variant={'psf', 'bg'} are equivalent, so both will have
                # variant={'bg'}.
                if variant in variants_done:
                    continue
                variants_done.add(variant)
                variant_dirname = '_'.join(sorted(variant))
                if len(variant) > 0:
                    variant_dirname = f'_{variant_dirname}'
                grouped_files = self.group_stage2_files_for_stage3(paths_in)

                # Only need to include the tile in prodname if it is needed to avoid
                # filename collisions for observations which use mosaicing. Most
                # observations just have a single tile per datset, so we can just use
                # the standard prodname in this case.
                keys_with_tiles = set(grouped_files.keys())
                keys_without_tiles = set(
                    (dither, match_key) for dither, tile, match_key in keys_with_tiles
                )
                include_tile_in_prodname = len(keys_with_tiles) != len(
                    keys_without_tiles
                )

                for (dither, tile, match_key), paths_in in grouped_files.items():
                    dirname = 'combined' if dither is None else f'd{dither}'
                    dirname = dirname + variant_dirname
                    output_dir = os.path.join(root_path, dir_out, dirname)
                    check_path(output_dir)

                    match_key_str = '_'.join(str(k) for k in match_key)
                    asn_path = os.path.join(
                        root_path,
                        dir_in,
                        f'{match_key_str}_dither-{dirname}_{tile}_asn.json',
                    )
                    prodname = 'Level3' + (
                        f'_{tile}' if include_tile_in_prodname else ''
                    )
                    self.write_asn_for_stage3(paths_in, asn_path, prodname=prodname)
                    paths_list.append((asn_path, output_dir))

            args_list = [
                (p, output_dir, kwargs) for p, output_dir in sorted(paths_list)
            ]
            runmany(
                self.reduction_spec3_fn,
                args_list,
                desc='stage3',
                **self.reduction_parallel_kwargs,
            )

    def reduction_spec3_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        asn_path, output_dir, kwargs = args
        Spec3Pipeline.call(
            asn_path,
            output_dir=output_dir,
            save_results=True,
            **kwargs,
        )

    def get_stage3_variant_paths_in(
        self, root_path: RootPath, variant: frozenset[str]
    ) -> tuple[frozenset[str], list[str]]:
        """
        Get list of input paths for a given variant for stage3.

        The returned variant may be a subset of the input variant if the variant
        contains a component produced after stage3. E.g. for MIRI, the PSF variants
        occur after stage3, so {'bg'} and {'bg', 'psf'} are equivalent when used as
        input to stage3.
        """
        dir_in, dir_out = self.step_directories['stage3']
        if variant == frozenset():
            paths = self.get_paths(root_path, dir_in, '*cal.fits')
        elif variant == frozenset({'bg'}):
            paths = self.get_paths(root_path, dir_in, 'bg', '*cal.fits')
        else:
            raise ValueError(f'Unknown variant: {variant}')
        return variant, paths

    def group_stage2_files_for_stage3(
        self, paths_in: list[str]
    ) -> dict[tuple[int | None, str, tuple[str | None, ...]], list[str]]:
        out: dict[tuple[int | None, str, tuple[str | None, ...]], list[str]] = {}
        dither_options = set()
        for p_in in paths_in:
            with fits.open(p_in) as hdul:
                hdr = hdul[0].header  #  type: ignore
            dither = int(hdr['PATT_NUM'])
            dither_options.add(dither)
            tile = hdr['ACT_ID']
            match_key = self.get_file_match_key_for_stage3(hdr)

            # Do both separated (dither, ...) and combined (None, ....) dithers
            for k in [
                (dither, tile, match_key),
                (None, tile, match_key),
            ]:
                out.setdefault(k, []).append(p_in)

        if len(dither_options) == 1:
            # If there is only a single dither, we can skip dither combination
            out = {
                (dither, tile, match_key): paths_in
                for (dither, tile, match_key), paths_in in out.items()
                if dither is not None
            }
        return out

    def get_file_match_key_for_stage3(
        self, path: str | fits.Header
    ) -> tuple[str | None, ...]:
        """
        Get a key used to match background and science files.
        """
        return self.get_file_match_key(path, self.stage3_file_match_hdr_keys)

    def write_asn_for_stage3(
        self, files: list[str], asn_path: str, prodname: str, **kwargs
    ) -> str:
        asn = asn_from_list(files, rule=DMS_Level3_Base, product_name=prodname)
        if 'bg' in kwargs:
            for bgfile in kwargs['bg']:
                asn['products'][0]['members'].append(
                    {'expname': bgfile, 'exptype': 'background'}
                )
        _, serialized = asn.dump()
        with open(asn_path, 'w', encoding='utf-8') as outfile:
            outfile.write(serialized)
        return asn_path

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
        for group_root_path in self.group_root_paths:
            paths_in = self.get_paths(
                group_root_path, dir_in, '*', '*_nav.fits', filter_variants=True
            )
            for p_in in paths_in:
                # get file path relative to the root of the group root path
                relpath = os.path.relpath(p_in, group_root_path)
                relpath = self.replace_path_part(relpath, -3, dir_out, old=dir_in)
                p_out = os.path.join(self.root_path, relpath)
                paths_dict.setdefault(p_out, []).append(p_in)

        args_list = [
            (paths_in, p_out, kwargs) for p_out, paths_in in paths_dict.items()
        ]
        runmany(
            self.desaturate_fn, args_list, desc='desaturate', **self.parallel_kwargs
        )

    def desaturate_fn(self, args: tuple[list[str], str, dict[str, Any]]) -> None:
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

    def despike_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        p_in, p_out, kwargs = args
        despike_data.despike_cube(p_in, p_out, **kwargs)

    # plot
    def run_plot(self, kwargs: dict[str, Any]) -> None:
        _, dir_out = self.step_directories['plot']
        for root_path in self.iterate_group_root_paths():
            paths = []
            for stage_dir in self.stage_directories_to_plot:
                paths_in = self.get_paths(
                    root_path, stage_dir, '*', '*_nav.fits', filter_variants=True
                )
                for p_in in paths_in:
                    dirname = os.path.dirname(
                        self.replace_path_part(
                            p_in, -3, os.path.join(dir_out, stage_dir), old=stage_dir
                        )
                    )
                    filename = self.get_plot_filename_prefix(p_in) + os.path.basename(
                        self.replace_path_suffix(p_in, '.png', old='.fits')
                    )
                    p_out = os.path.join(dirname, filename)
                    paths.append((p_in, p_out))
            args_list = [(p_in, p_out, kwargs) for p_in, p_out in paths]
            runmany(self.plot_fn, args_list, desc='plot', **self.parallel_kwargs)

    def plot_fn(self, args: tuple[str, str, dict[str, Any]]) -> None:
        p_in, p_out, kwargs = args
        jwst_summary_plots.make_summary_plot(p_in, p_out, **kwargs)

    # pylint: disable-next=unused-argument
    def get_plot_filename_prefix(self, path: str) -> str:
        """
        Get filename prefix to use for plots/animations e.g. '1A_' for MIRI.
        """
        return ''

    # animate
    def run_animate(self, kwargs: dict[str, Any]) -> None:
        _, dir_out = self.step_directories['animate']
        paths_dict: dict[str, list[str]] = {}
        for stage_dir in self.stage_directories_to_plot:
            # use d* as dither path as we don't want to use combined dithers here
            paths_in = self.get_paths(
                self.root_path, stage_dir, 'd*', '*_nav.fits', filter_variants=True
            )
            for p_in in paths_in:
                variant_dir = pathlib.Path(p_in).parts[-2].split('_')[1:]
                variant_dir = 'dither_comparison_' + '_'.join(variant_dir)
                variant_dir = variant_dir.rstrip('_')
                dirname = os.path.join(self.root_path, dir_out, stage_dir, variant_dir)
                filename = self.get_plot_filename_prefix(p_in) + os.path.basename(
                    self.replace_path_suffix(p_in, '.mp4', old='.fits')
                )
                p_out = os.path.join(dirname, filename)
                paths_dict.setdefault(p_out, []).append(p_in)
        args_list = [
            (paths_in, p_out, kwargs) for p_out, paths_in in paths_dict.items()
        ]
        runmany(self.animate_fn, args_list, desc='animate', **self.parallel_kwargs)

    def animate_fn(self, args: tuple[list[str], str, dict[str, Any]]) -> None:
        paths_in, p_out, kwargs = args
        jwst_summary_animation.make_animation(paths_in, p_out, **kwargs)


def get_pipeline_argument_parser(
    pipeline_class: Type[Pipeline],
    step_descriptions: str,
    cli_examples: str,
) -> argparse.ArgumentParser:
    """
    Get a CLI argument parser for a given pipeline class.

    The returned parser instance can be further customised before calling
    `run_pipeline(**vars(parser.parse_args()))` to run the pipeline.
    """
    suffix = step_descriptions + '\n' + cli_examples
    name = pipeline_class.get_instrument()
    parser = argparse.ArgumentParser(
        description=(
            f'Full JWST {name} IFU pipeline including the standard reduction from '
            'stage0 to stage3, and custom pipeline steps for additional cleaning and '
            'data visualisation. For more customisation, this script can be imported '
            'and run in Python (see the source code for mode details).\n\n'
            'The following steps are run in the full pipeline:' + suffix
        ),
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
        '--groups-to-use',
        type=str,
        help="""Comma-separated list of groups to keep. For example, `1,2,3,4` will
            keep the first four groups. If unspecified, all groups will be kept.""",
    )
    parser.add_argument(
        '--background_subtract',
        '--background-subtract',
        action=argparse.BooleanOptionalAction,
        default='both',
        help="""Toggle background subtractio of the data. If unspecified, then versions
        with and without background subtraction will be created. Background subtraction
        requires --background_path to be specified.""",
    )
    parser.add_argument(
        '--background_path',
        '--background-path',
        type=str,
        help="""Path to directory containing background data. For example, if
            your `root_path` is `/data/uranus/lon1`, the `background_path` may
            be `/data/uranus/background`. Note that background subtraction will
            require the background data to be already reduced to `stage1`. If no 
            `background_path` is specified (the default), then no background subtraction
            will be performed.""",
    )
    parser.add_argument(
        '--basic_navigation',
        '--basic-navigation',
        action='store_true',
        help="""Use basic navigation, and only save RA and Dec backplanes (e.g. useful
            for small bodies). By default, full navigation is performed, generating a
            full set of coordinate backplanes (lon/lat, illumination angles etc.). Using
            basic navigation automatically skips the animation step. This
            is mainly useful if you get SPICE errors when navigating the data.""",
    )
    parser.add_argument(
        '--step_kwargs',
        '--step-kwargs',
        '--kwargs',
        type=str,
        help="""JSON string containing keyword arguments to pass to individual pipeline
            steps. For example, 
            `--kwargs '{"stage3": {"steps": {"outlier_detection": {"snr": "30.0 24.0", "scale": "1.3 0.7"}}}, "plot": {"plot_brightest_spectrum": true}}'` 
            will pass the custom arguments to the stage3 and plot steps.
            """,
    )
    parser.add_argument(
        '--skip_steps',
        '--skip-steps',
        nargs='+',
        type=str,
        help="""List of steps to skip. This is generally only useful if you are
            re-running part of the pipeline. Multiple steps can be passed as a
            space-separated list. For example, `--skip_steps flat despike`.""",
    )
    parser.add_argument(
        '--start_step',
        '--start-step',
        type=str,
        help="""Convenience argument to add all steps before `start_step` to 
            `skip_steps`.""",
    )
    parser.add_argument(
        '--end_step',
        '--end-step',
        type=str,
        help="""Convenience argument to add all steps steps after `end_step` to 
            `skip_steps`.""",
    )
    return parser
