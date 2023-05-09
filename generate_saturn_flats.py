#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate flat fields for Saturn observations (Fletcher et al. 2023).
"""
import os
import tqdm
import construct_flat_field
import tools

TILES = ['SATURN-15N', 'SATURN-45N']
CHANNEL_NUMBERS = ['1', '2', '3', '4']
CHANNEL_BANDS = ['short', 'medium', 'long']
FRINGES = ['_fringe', '']
DITHERS = ['1', '2', '3', '4']

ROOT_PATH = '/data/nemesis/jwst/MIRI_IFU/Saturn_2022nov13'
PATTERN_IN = os.path.join(
    ROOT_PATH,
    '{dataset}',
    'stage3_desaturated',
    'd{dither}{fringe}_nav',
    'Level3_ch{channel_number}-{channel_band}_s3d_nav.fits',
)
PATTERN_OUT = os.path.join(
    ROOT_PATH,
    'flat_field',
    '{dataset}',
    'ch{channel_number}-{channel_band}_constructed_flat{fringe}.fits',
)

LAT_BIN_SIZE_FACTOR = 0.5
BIN_ASPECT = 4


def main():
    for kwargs in tqdm.tqdm(
        tools.all_combinations(
            fringe=FRINGES, channel_number=CHANNEL_NUMBERS, channel_band=CHANNEL_BANDS
        ),
        desc='Generating flats',
    ):
        do_channel(**kwargs)


def do_channel(channel_number: str, channel_band: str, fringe: str):
    try:
        for tile in tqdm.tqdm(
            TILES,
            desc=f'Creating {channel_number}-{channel_band}{fringe} flats',
            leave=False,
        ):
            do_tile(channel_number, channel_band, fringe, tile)
        merge_flats(channel_number, channel_band, fringe)
    except FileNotFoundError:
        print(f'No data found for ch{channel_number}-{channel_band}{fringe}')


def do_tile(channel_number: str, channel_band: str, fringe: str, dataset: str):
    kw = dict(
        channel_number=channel_number,
        channel_band=channel_band,
        dataset=dataset,
        fringe=fringe,
    )
    paths = [PATTERN_IN.format(**kw, dither=dither) for dither in DITHERS]
    p_out = PATTERN_OUT.format(**kw)
    construct_flat_field.do_tile(
        paths,
        p_out,
        channel=channel_number,
        band=channel_band,
        fringe=fringe,
        dataset=dataset,
        lat_bin_size_factor=LAT_BIN_SIZE_FACTOR,
        bin_aspect=BIN_ASPECT,
    )


def merge_flats(channel_number: str, channel_band: str, fringe: str):
    kw = dict(
        channel_number=channel_number,
        channel_band=channel_band,
        fringe=fringe,
    )
    paths = [PATTERN_OUT.format(**kw, dataset=ds) for ds in TILES]
    p_out = PATTERN_OUT.format(**kw, dataset='merged')
    construct_flat_field.merge_flats(paths, p_out)


if __name__ == '__main__':
    main()
