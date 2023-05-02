"""
File path record.

Central location to modify and save paths in file system.
"""
import os
import sys


def generate_path(*args: str) -> str:
    """
    Create consistently formatted file/directory path string.

    Parameters
    ----------
    args : str
        Series of strings representing parts of path

    Returns
    -------
    str
    """
    path = os.path.sep.join(args)
    path = os.path.normpath(path)
    return path


_ALICE = sys.platform == 'linux'

if _ALICE:
    DATA_ROOT = generate_path('/data/nemesis/jwst')
else:
    DATA_ROOT = generate_path('/Users/ortk1/Dropbox/science/jwst_data')


def jwst(*args: str) -> str:
    return generate_path(DATA_ROOT, *args)


def saturn(*args: str) -> str:
    return jwst('MIRI_IFU/Saturn_2022nov13', *args)
