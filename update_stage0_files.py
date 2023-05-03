#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import tools
import os
import tqdm
import glob
import urllib.parse
import urllib.request
from astropy.io import fits
import sys

MAST_URL = (
    'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/{filename}'
)


def main(*args: str) -> None:
    for directory in args:
        update_directory(directory)


def update_directory(directory_path: str, **kw) -> None:
    """
    Update all files in a given directory.

    Args:
        directory_path: Path to directory.
    """
    files = sorted(glob.glob(os.path.join(directory_path, '*.fits')))
    print(f'Updating {len(files)} files in {directory_path}...')
    for file in files:
        update_file(file, **kw)
    print(f'Successfully updated {len(files)} files in {directory_path}.')


def update_file(path: str, remove_old: bool = True) -> None:
    """
    Downloads new version of stage0 file from MAST.

    Args:
        path: Path of local stage0 file to update.
        remove_old: Remove the old file if the filename has changed.
    """
    with fits.open(path) as hdul:
        primary_header: fits.Header = hdul['PRIMARY'].header  # type: ignore
        filename = str(primary_header['FILENAME'])
    url = MAST_URL.format(filename=filename)

    path_out = os.path.join(os.path.split(path)[0], filename)
    download_file(url, path_out)

    if remove_old and path != path_out:
        print(f'Removed old file {os.path.basename(path)}')
        os.remove(path)


def download_file(url: str, local_path: str) -> None:
    """
    Download file to local system.

    Args:
        url: URL of file file.
        local_path: File path on local system.
    """
    tools.check_path(local_path)

    # download to temp file so don't get issues from partial downloads being killed
    temp_path = local_path + '.temp'
    urllib.request.urlretrieve(
        url, temp_path, reporthook=DownloadProgressBar(os.path.basename(local_path))
    )

    # once fully downloaded, we can safely move the temp file to the desired path
    os.replace(temp_path, local_path)


class DownloadProgressBar:
    """
    Shows download progress with tqdm
    """

    def __init__(self, filename: str):
        self.pbar = None
        self.previous_downloaded = 0
        self.filename = filename

    def __call__(self, block_num, block_size, total_size):
        if not self.pbar:
            self.pbar = tqdm.tqdm(
                total=total_size,
                unit_scale=True,
                unit='B',
                unit_divisor=1024,
                desc=self.filename,
            )
        downloaded = block_num * block_size
        change = downloaded - self.previous_downloaded
        self.previous_downloaded = downloaded
        self.pbar.update(change)


if __name__ == '__main__':
    main(*sys.argv[1:])
