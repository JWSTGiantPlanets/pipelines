#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JWST reduction script for MIRI MRS data - uses standard pipeline to reduce to stage3.

Runs stages 1-3 of the JWST pipeline with optional defringing step after stage 2. Does
not combine dithers. Based on Jake's jwstpipeline_singledither_v2.py.

Before running this script, you likely want to set the following environment variables
(e.g. in your ~/.bash_profile script) ::

    export CRDS_PATH="path/to/crds_cache"
    export CRDS_SERVER_URL="https://jwst-crds.stsci.edu"

This pipeline can be run from the command line (see CLI_EXAMPLES below) or you can 
import this module and use the JWSTReduction class directly. ::

    import reduce_jwst
    pipeline = reduce_jwst.JWSTReduction('data/saturn/SATURN-15N')
    pipeline.run(defringe=True)
"""
CLI_EXAMPLES = """examples:

# Print a help message
python3 reduce_jwst.py -h

# Run all stages of the pipeline on the data in data/saturn/SATURN-15N/stage0
python3 reduce_jwst.py data/saturn/SATURN-15N

# Run stages 2 and 3 of the pipeline
python3 reduce_jwst.py data/saturn/SATURN-15N -s 2 3

# Run stages 2 and 3 of the pipeline and defringes the data
python3 reduce_jwst.py data/saturn/SATURN-15N -s 2 3 -d

# Run all stages of the pipeline in parallel using half of the available cores
python3 reduce_jwst.py data/saturn/SATURN-15N -p 0.5
"""
import argparse
import datetime
import functools
import glob
import multiprocessing
import os
from typing import Any, Callable

import numpy as np
from astropy.io import fits

import parallel_tools
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.pipeline import Detector1Pipeline, Spec2Pipeline, Spec3Pipeline
from jwst.residual_fringe import ResidualFringeStep


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='JWST standard reduction pipeline (stage0 to stage3)',
        epilog=CLI_EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        argument_default=argparse.SUPPRESS,
    )
    parser.add_argument(
        'root_path',
        type=str,
        help='Path of location of stage0 files',
    )
    parser.add_argument(
        '--stages',
        '-s',
        nargs='*',
        type=float,
        help='Stages to run',
    )
    parser.add_argument(
        '--defringe',
        action='store_true',
        help='Run defringe steps',
    )
    parser.add_argument(
        '--parallel',
        nargs='?',
        const=1,
        type=float,
        help='Fraction of cores to use when multiprocessing, set to 0 to run serially or 1 to use all cores',
    )
    parser.add_argument(
        '--parallel_timeout',
        type=float,
        help='Timeout in seconds for parallel tasks, set to 0 to disable timeout',
    )
    args = vars(parser.parse_args())
    JWSTReduction(**args, _run_immediately=True)


class JWSTReduction:
    FILENAME_PREFIX = ''
    STAGE0_SUFFIX = ''

    def __init__(
        self,
        root_path: str,
        parallel: float | bool = False,
        parallel_timeout: float | None = 60 * 60,
        **kwargs,
    ) -> None:
        """
        JWST reduction pipeline.

        Example usage to run stages 2 and 3 with defringing on all cores: ::

            import reduce_jwst
            pipeline = reduce_jwst.JWSTReduction('data/saturn/SATURN-15N')
            pipeline.run([2,3], defringe=True, parallel=1)

        Args:
            root_path: Path to directory containing stage0 files.
            parallel: Fraction of CPU cores to use when multiprocessing. Set to 0 to run
                serially, 1 to use all cores, 0.5 to use half of the cores, etc.
            parallel_timeout: Timeout in seconds for parallel tasks. This is the
                maximum average time for a batch of parallel tasks to complete. Set to
                None or 0 to disable timeout. Timed out jobs will be reattempted
                serially.

        Raises:
            FileNotFoundError: If no stage0 directory is found.
        """
        parallel = float(parallel)
        root_path = os.path.normpath(root_path)
        self.root_path = root_path
        if parallel_timeout == 0:
            parallel_timeout = None
        self.parallel_timeout = parallel_timeout
        self.parallel_start_delay = 10

        self.stage0_dir = os.path.join(root_path, 'stage0')
        self.stage1_dir = os.path.join(root_path, 'stage1')
        self.stage2_dir = os.path.join(root_path, 'stage2')
        self.stage3_dir = os.path.join(root_path, 'stage3')

        if not os.path.exists(self.stage0_dir):
            raise FileNotFoundError(
                f'stage0 directory {self.stage0_dir!r} does not exist'
            )

        self._max_processors = parallel_tools.get_max_processors(parallel)

        if '_run_immediately' in kwargs:
            kwargs.pop('_run_immediately')
            self.run(**kwargs)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}({self.root_path!r})'

    def run(self, stages: list | None = None, defringe: bool = False) -> None:
        """
        Run pipeline steps.

        Args:
            stages: List of stages to run. Defaults to `[1, 2, 3]` to run all stages.
            defringe: Toggle running the defringe step and using defringed data in
                stage 3.
        """
        if stages is None:
            stages = [1, 2, 3]
        stages = [(s.replace('stage', '') if isinstance(s, str) else s) for s in stages]
        stages = [float(s) for s in stages]

        log('Running pipeline stages:', ', '.join(format(s, 'g') for s in stages))
        log('Defringe:', defringe, time=False)
        log('Path:', self.root_path, time=False)
        log(
            'Running',
            'serially' if self._max_processors == 1 else 'in parallel',
            f'on {self._max_processors}/{multiprocessing.cpu_count()} cores',
            time=False,
        )
        print()

        if 1 in stages:
            self.run_stage1()

        if 2 in stages:
            self.run_stage2()

        if 2 in stages or 2.5 in stages:
            if defringe:
                self.run_defringe()

        if 3 in stages:
            self.run_stage3(defringe)

    def run_stage1(self) -> None:
        log('Running stage 1 pipeline')
        raw_files = sorted(
            glob.glob(
                os.path.join(
                    self.stage0_dir, f'{self.FILENAME_PREFIX}*{self.STAGE0_SUFFIX}.fits'
                )
            )
        )
        if not os.path.exists(self.stage1_dir):
            os.makedirs(self.stage1_dir)
        self._runmany(self._rundet1, raw_files)
        log('Finished stage 1 pipeline\n')

    def run_stage2(self) -> None:
        log('Running stage 2 pipeline')
        rate_files = sorted(
            glob.glob(
                os.path.join(self.stage1_dir, f'{self.FILENAME_PREFIX}*rate.fits')
            )
        )
        for path in rate_files:
            with fits.open(path) as hdul:
                hdul = fits.open(path)
                hdul[1].header['SRCTYPE'] = 'EXTENDED'  # type: ignore
                hdul.writeto(path, overwrite=True)
        if not os.path.exists(self.stage2_dir):
            os.makedirs(self.stage2_dir)
        self._runmany(self._runspec2, rate_files)
        log('Finished stage 2 pipeline\n')

    def run_defringe(self) -> None:
        log('Running residual fringe step')
        cal_files = sorted(
            glob.glob(os.path.join(self.stage2_dir, f'{self.FILENAME_PREFIX}*cal.fits'))
        )
        self._runmany(self._runfringe, cal_files)
        print('Finished residual fringe step\n')

    def run_stage3(self, defringe: bool) -> None:
        log('Running stage 3 pipeline')
        filename_suffix = 'residual_fringe.fits' if defringe else 'cal.fits'
        calfiles = sorted(
            glob.glob(
                os.path.join(
                    self.stage2_dir, f'{self.FILENAME_PREFIX}*{filename_suffix}'
                )
            )
        )

        sorted_files = self._sort_calfiles(calfiles)
        for dither, dither_files in sorted_files.items():
            log(f'Processing dither {dither}')

            output_path = os.path.join(self.stage3_dir, f'd{dither}')
            if defringe:
                output_path = output_path + '_fringe'
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            asnlist = []
            for name, files in dither_files.items():
                if len(files) > 0:
                    asn_path = os.path.join(self.stage2_dir, f'l3asn-{name}.json')
                    asnlist.append(asn_path)
                    self._writel3asn(files, asn_path, 'Level3')

            fn = functools.partial(self._runspec3, output_path=output_path)
            self._runmany(fn, asnlist)

        log('Finished stage 3 pipeline\n')

    def _rundet1(self, path: str) -> None:
        with fits.open(path) as hdul:
            ngroups = hdul['PRIMARY'].header['NGROUPS']  # type: ignore

        det1 = Detector1Pipeline()
        det1.output_dir = self.stage1_dir  # type: ignore
        det1.dq_init.skip = False  # type: ignore
        det1.saturation.skip = False  # type: ignore
        det1.superbias.skip = False  # type: ignore
        det1.refpix.skip = False  # type: ignore
        det1.linearity.skip = False  # type: ignore
        det1.persistence.skip = False  # type: ignore
        det1.dark_current.skip = False  # type: ignore

        if ngroups <= 3:
            det1.firstframe.skip = True  # type: ignore
            det1.lastframe.skip = True  # type: ignore
            det1.rscd.skip = True  # type: ignore
            det1.jump.skip = True  #  type: ignore

        det1.save_results = True
        det1(path)

    def _runspec2(self, path: str) -> None:
        spec2 = Spec2Pipeline()
        spec2.output_dir = self.stage2_dir  # type: ignore
        spec2.assign_wcs.skip = False  # type: ignore
        spec2.bkg_subtract.skip = True  # type: ignore
        spec2.flat_field.skip = False  # type: ignore
        spec2.srctype.skip = False  # type: ignore
        spec2.straylight.skip = False  # type: ignore
        spec2.fringe.skip = False  # type: ignore
        spec2.photom.skip = False  # type: ignore
        spec2.cube_build.skip = True  # type: ignore
        spec2.extract_1d.skip = True  # type: ignore
        spec2.save_results = True
        spec2(path)

    def _runfringe(self, path: str) -> None:
        rf1 = ResidualFringeStep()
        rf1.skip = False
        rf1.save_results = True
        rf1.output_dir = self.stage2_dir  # type: ignore
        rf1(path)

    def _runspec3(self, path: str, output_path: str) -> None:
        crds_config = Spec3Pipeline.get_config_from_reference(path)
        spec3 = Spec3Pipeline.from_config_section(crds_config)
        spec3.output_dir = output_path  # type: ignore
        spec3.save_results = True  # type: ignore
        spec3.assign_mtwcs.skip = False  # type: ignore
        spec3.master_background.skip = True  # type: ignore
        spec3.outlier_detection.skip = False  # type: ignore
        spec3.mrs_imatch.skip = False  # type: ignore
        spec3.cube_build.skip = False  # type: ignore
        spec3.extract_1d.skip = True  # type: ignore
        spec3.cube_build.output_type = 'band'  # type: ignore
        spec3.cube_build.coord_system = 'ifualign'  # type: ignore
        spec3(path)

    @staticmethod
    def _sort_calfiles(paths: list[str]) -> dict[int, dict[str, list[str]]]:
        channel = []
        band = []
        dither = []
        for p in paths:
            with fits.open(p) as hdul:
                hdr = hdul[0].header  # type: ignore
            channel.append(hdr['CHANNEL'])
            band.append(hdr['BAND'])
            dither.append(hdr['PATT_NUM'])
        channel = np.array(channel)
        band = np.array(band)
        dither = np.array(dither)
        output = {}
        path_array = np.array(paths)
        dithers = sorted(set(dither))
        for d in dithers:
            output[d] = {}
            for c in ['12', '34']:
                for b, abc in zip(['SHORT', 'MEDIUM', 'LONG'], ['A', 'B', 'C']):
                    indx = np.where((channel == c) & (band == b) & (dither == d))
                    output[d][c + abc] = list(path_array[indx])
        return output

    @staticmethod
    def _writel3asn(files, asnfile, prodname, **kwargs):
        asn = afl.asn_from_list(files, rule=DMS_Level3_Base, product_name=prodname)
        if 'bg' in kwargs:
            for bgfile in kwargs['bg']:
                asn['products'][0]['members'].append(
                    {'expname': bgfile, 'exptype': 'background'}
                )
        _, serialized = asn.dump()
        with open(asnfile, 'w') as outfile:
            outfile.write(serialized)

    def _runmany(self, step: Callable[[str], Any], filepaths: list[str]) -> None:
        parallel_tools.runmany(
            step,
            filepaths,
            num_processors=self._max_processors,
            timeout=self.parallel_timeout,
            start_delay=self.parallel_start_delay,
            parallel_job_kw=dict(
                caught_error_wait_time=60,
                caught_error_wait_time_frac=1,
                caught_error_wait_time_max=600,
            ),
        )


def log(*messages: Any, time: bool = True) -> None:
    prefix = datetime.datetime.now().strftime('%H:%M:%S') if time else ' ' * 8
    print(prefix, *messages, flush=True)


if __name__ == '__main__':
    main()
