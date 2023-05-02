import datetime
import math
import multiprocessing
import multiprocessing.pool
import os
import time
from typing import Any, Callable, TypeVar

import tqdm

from tools import log

T = TypeVar('T')


def runmany(
    function: Callable[[T], Any],
    args_list: list[T],
    *,
    num_processors: int | None = None,
    parallel_frac: float | bool = False,
    timeout: float | None = None,
    start_delay: float | None = None,
    retry_parallel_errors: bool = True,
    parallel_job_kw: dict[str, Any] | None = None,
    **tqdm_kw,
) -> None:
    """
    Function to run many jobs in serial/parallel depending on `parallel` argument.

    Args:
        function: Function to run on each argument.
        args_list: List of arguments to pass to `function`.
        num_processors: Maximum number of processors to use. If `None`, will
            calculate based on `parallel_frac`.
        parallel_frac: Fraction of processors to use if `num_processors` is `None`. Set
            to `False` to run serially (i.e. completely avoid using parallel
            processing). Note that if running serially, then the `timeout` and 
            `start_delay` arguments are ignored.
        timeout: Average timeout for each job in seconds. If `None`, will not timeout.
            This timeout is performed on the entire batch of jobs, so long jobs run at
            the start of the batch may cause later jobs to timeout if the entire queue
            exceeds the timeout. The total timeout for the entire batch is calculated as
            `timeout * math.ceil(len(args_list) / num_processors)`. If
            `retry_parallel_errors` is `True`, then timed out jobs will be rerun
            serially.
        start_delay: Delay between starting each job in seconds. This can be used to
            help avioid any race conditions that may occur when starting multiple jobs
            at the same time. Set to `None` to have no delay.
        retry_parallel_errors: If `True`, will retry jobs that fail in parallel by
            running serially once all parallel jobs have completed (or failed).
        parallel_job_kw: Keyword arguments to pass to `parallel_job`.
        tqdm_kw: Keyword arguments to pass to tqdm progress bar when running serially.
    """
    if num_processors is None:
        num_processors = get_max_processors(parallel_frac)

    if num_processors == 1:
        for a in tqdm.tqdm(args_list, **tqdm_kw):
            function(a)
    else:
        log(
            f'Processing {len(args_list)} jobs in parallel using {num_processors} cores...'
        )
        if timeout is not None:
            timeout *= max(1, math.ceil(len(args_list) / num_processors))
        dtm_start = datetime.datetime.now()
        parallel_args = make_parallel_args(function, args_list)
        with multiprocessing.get_context('spawn').Pool(num_processors) as p:
            jobs: list[multiprocessing.pool.ApplyResult[bool]] = []
            for a in parallel_args:
                jobs.append(
                    p.apply_async(
                        parallel_job,
                        (a,),
                        dict(
                            catch_errors=retry_parallel_errors,
                            **(parallel_job_kw or {}),
                        ),
                    )
                )
                if start_delay:
                    time.sleep(start_delay)
            
            to_retry = wait_for_jobs(
                jobs, timeout=timeout, catch_errors=retry_parallel_errors
            )

        if any(to_retry):
            args_list_retry = [a for a, r in zip(parallel_args, to_retry) if r]
            log(f'Rerunning {len(args_list_retry)} failed jobs serially...')
            for a in args_list_retry:
                parallel_job(a, catch_errors=False)
            log(f'Completed {len(args_list_retry)} rerunning jobs serially')
        
        dtm_end = datetime.datetime.now()
        log(f'Completed {len(args_list)} jobs in {dtm_end - dtm_start}')


def get_max_processors(parallel_frac: float = 1) -> int:
    try:
        # Get the number of processors per node when running on ALICE
        num_cores = int(os.environ['PBS_NUM_PPN'])
    except (KeyError, ValueError):
        num_cores = multiprocessing.cpu_count()
    return max(int(num_cores * parallel_frac), 1)


def wait_for_jobs(
    jobs: list[multiprocessing.pool.ApplyResult],
    timeout: float | None,
    catch_errors: bool,
) -> list[bool]:
    results: list[bool] = []
    t0 = time.time()
    for idx, job in enumerate(jobs):
        try:
            if timeout is not None:
                # Reduce timeout by elapsed time
                timeout = max(0, timeout - (time.time() - t0))
                t0 = time.time()
            res = job.get(timeout=timeout)
        except multiprocessing.TimeoutError:
            if catch_errors:
                log(
                    f'>> Job {idx+1}/{len(jobs)} timed out - will retry serially later...'
                )
                res = True
            else:
                log(f'>> Job {idx+1}/{len(jobs)} timed out')
                raise
        results.append(res)
    return results


def make_parallel_args(
    fn: Callable[[T], Any], args_list: list[T]
) -> list[tuple[Callable[[T], Any], T, int, int]]:
    """
    Helper function for `runmany` to create arguments for parallel jobs.
    """
    return [(fn, args, idx, len(args_list)) for idx, args in enumerate(args_list)]


def parallel_job(
    job_args: tuple[Callable[[T], Any], T, int, int],
    catch_errors: bool,
    caught_error_wait_time: float = 0,
    caught_error_wait_time_frac: float = 0,
    caught_error_wait_time_max: float | None = None,
) -> bool:
    """
    Helper function for `runmany` to run a parallel job.

    Args:
        job_args: Tuple of function, arguments, job index, and total number of jobs.
        catch_errors: If `True`, will catch errors and return `True` to enable retrying
            failed jobs serially.
        caught_error_wait_time: Time to wait after catching an error before returning.
            This can be used to help prevent errors building up by multiple tasks
            sequentially failing due to race conditions.
        caught_error_wait_time_frac: Fraction of time taken by job to wait after
            catching an error before returning.
        caught_error_wait_time_max: Maximum time to wait after catching an error before
            returning.
    """
    fn, args, idx, total = job_args
    log(f'>> Job {idx+1}/{total} starting...')
    d0 = datetime.datetime.now()
    try:
        fn(args)
    except Exception as e:
        if catch_errors:
            log(f'>> Job {idx+1}/{total} failed - will retry serially later...')
            log(f'>> Caught error: {e!r}', time=False)
            d1 = datetime.datetime.now()
            t = (
                caught_error_wait_time
                + caught_error_wait_time_frac * (d1 - d0).total_seconds()
            )
            if caught_error_wait_time_max is not None:
                t = min(t, caught_error_wait_time_max)
            if t > 0:
                log(
                    f'>> Waiting {t}s before next job...',
                    time=False,
                )
                time.sleep(t)
            return True
        log(f'>> Job {idx+1}/{total} failed with arguments: {args}')
        raise
    d1 = datetime.datetime.now()
    log(f'>> Job {idx+1}/{total} completed in {d1-d0}')
    return False
