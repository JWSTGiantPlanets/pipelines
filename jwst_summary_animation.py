#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create useful summary animations of JWST NIRSPEC and MIRI observations.

The animations scroll through the wavelengths of the provided cubes, showing multiple
different projections of the data. The animations can compare multiple different cubes,
so can be used to e.g. compare all 4 dithers of a single observation. Each animation 
takes ~2 minutes to generate.

The input files to this script must be navigated using the JWSTSolarSystemPointing code:
https://github.com/JWSTGiantPlanets/NIRSPEC-Toolkit/blob/main/navigate_jwst_observations.py

See also jwst_summary_plots.py

This can be run from the command line (see CLI_EXAMPLES below) or from Python (see
docstring within each `make_animation()` below).
"""
CLI_EXAMPLES = """examples:

python3 jwst_summary_animation.py -h
# Prints a help message

python3 jwst_summary_animation.py obs1/d?_nav/*ch1-medium*.fits animation.mp4
# Creates animation of all dithers of channel 1B off the obs1 tile

python3 jwst_summary_animation.py --wavelength_step 1 --fps 30 obs1/d?_nav/*ch1-medium*.fits animation.mp4 
# Creates animation where each frame is a single wavelength at 30 FPS

python3 jwst_summary_animation.py --radius_factor 3 obs1/d?_nav/*ch1-medium*.fits animation.mp4 
# Creates animation where the target takes up 1/3 of the context plot (e.g. to see rings)
"""
import argparse
import os
import sys
from typing import Literal

import matplotlib.animation as mpanim
import matplotlib.pyplot as plt
import matplotlib.style as mstyle
import numpy as np
import planetmapper
import tqdm
from astropy.io import fits
from matplotlib.cm import ScalarMappable

import tools
import zonal_average

PLANETMAPPER_KW = {}
if sys.platform == 'linux' and 'PLANETMAPPER_KERNEL_PATH' not in os.environ:
    # Try to set the kernel path automatically when running on ALICE
    PLANETMAPPER_KW['kernel_path'] = '/data/nemesis/jwst/scripts/kernels'

OBS_CMAP = 'viridis'
TB_CMAP = 'inferno'

OBS_CMAP_DIFF = 'BrBG_r'
TB_CMAP_DIFF = 'RdBu_r'

FIGURE_BACKGROUND_COLOR = 'k'
AXES_BACKGROUND_COLOR = '0.075'
WIREFRAME_COLOR = '0.5'


def make_animation(
    input_paths: list[str],
    output_path: str,
    fps: int = 20,
    wavelength_step: int = 5,
    map_degree_interval: float = 1,
    map_interpolation: Literal['nearest', 'linear', 'quadratic', 'cubic'] = 'linear',
    clip_percentile: float = 0.2,
    clip_allowed_range: float = 0.2,
    radius_factor: float = 1.25,
    dpi: int = 150,
    title_suffix: str | list[str] | None = None,
    zonal_average_path: str | None = None,
    progress_bar: bool = True,
    print_output: bool = True,
) -> None:
    """
    Create animation of observed cubes. Example usage ::

        import jwst_summary_animation
        import glob

        input_paths = sorted(glob.glob('obs1/d?_nav/Level3_ch1-medium_s3d_nav.fits'))
        jwst_summary_animation.make_animation(input_paths, 'animations/ch1-medium.mp4')

    Args:
        input_paths: List of paths to navigated FITS cubes. This should generally be a
            list of 4 paths corresponding to the 4 dithers within a single observation.
        output_path: Path of output animation
        fps: Frames per second of output animation.
        wavelength_step: Number of wavelengths averaged per frame of animation.
        map_degree_interval: Degree interval of plotted map.
        map_interpolation: Interpolation used when projecting mapped data. This can be
             any of `'nearest'`, `'linear'`, `'quadratic'` or `'cubic'`.
        clip_percentile: Percentile of to set colour limits (to identify outliers). See
            `get_limits()` for where this is used.
        clip_allowed_range: Tolerance factor for identifying outliers outside the
            percentile range.  See `get_limits()` for where this is used.
        radius_factor: Size of context plot limits in multiples of the target's
            equatorial diameter.
        dpi: Resolution of output animation. The figure size is 12x12 inches, so the
            default DPI of 150 creates a 1800x1800 pixel animation.
        title_suffix: Additional suffix to add to the title of the animation.
        zonal_average_path: Path of zonal average data to plot deviations. If this is
            `None` (the default), then the absolute values of the data are plotted as
            normal.
        progress_bar: Whether to show a progress bar.
        print_output: Whether to print progress info to the console.
    """
    if print_output:
        print(
            'Loading and preparing data for {}...'.format(os.path.split(output_path)[1])
        )
    observations: list[planetmapper.Observation] = []
    headers: list[tuple[fits.Header, fits.Header]] = []
    data_wavelengths: list[np.ndarray] = []
    lat_images: list[np.ndarray] = []

    title = None
    for p in input_paths:
        with fits.open(p) as hdul:
            primary_header: fits.Header = hdul['PRIMARY'].header  # type: ignore
            data_header: fits.Header = hdul['SCI'].header  # type: ignore
            ra_img = hdul['RA'].data  # type: ignore
            dec_img = hdul['DEC'].data  # type: ignore
            lat_img = hdul['LAT_GRAPHIC'].data  # type: ignore
            wl = tools.get_wavelengths(data_header)

        obs = planetmapper.Observation(p, **PLANETMAPPER_KW)
        obs.disc_from_wcs(suppress_warnings=True)

        # Apply any RA/Dec offsets applied to the data when navigating
        dra = np.nanmean(ra_img - obs.get_ra_img()) * 60 * 60
        ddec = np.nanmean(dec_img - obs.get_dec_img()) * 60 * 60
        obs.add_arcsec_offset(dra, ddec)

        observations.append(obs)
        headers.append((primary_header, data_header))
        data_wavelengths.append(wl)
        lat_images.append(lat_img)

        instrument = primary_header['INSTRUME']
        title_parts = [
            primary_header['TARGPROP'],
            primary_header['INSTRUME'],
        ]
        if instrument == 'MIRI':
            title_parts.append(
                'Channel {} {}'.format(
                    primary_header['CHANNEL'], primary_header['BAND']
                )
            )
        if instrument == 'NIRSPEC':
            title_parts.append(
                '{} {} {}'.format(
                    primary_header['DETECTOR'],
                    primary_header['FILTER'],
                    primary_header['GRATING'],
                )
            )
        if zonal_average_path is not None:
            title_parts.append('Deviation from zonal average')
        if title_suffix is not None:
            if isinstance(title_suffix, str):
                title_suffix = [title_suffix]
            title_parts.extend(title_suffix)
        title_loc = '  |  '.join(title_parts)  #  type: ignore
        if title is None:
            title = title_loc
        else:
            assert title_loc == title
    assert title is not None
    for w in data_wavelengths:
        assert np.allclose(data_wavelengths[0], w, equal_nan=True)
    wavelengths = data_wavelengths[0]

    with tools.ignore_warnings(
        'invalid value encountered in log', 'divide by zero encountered'
    ):
        tb_cubes = [
            tools.brightness_temperature(obs.data, data_header)
            for obs, (primary_header, data_header) in zip(observations, headers)
        ]
    obs_cubes = [obs.data for obs in observations]
    mask_images = [np.isfinite(l) for l in lat_images]

    spectra = [np.nanmean(tb_cube, axis=(1, 2)) for tb_cube in tb_cubes]
    spectrum = np.nanmedian(spectra, axis=0)

    if zonal_average_path is not None:
        cube_zonal, hdr_zonal = fits.getdata(zonal_average_path, header=True)  # type: ignore
        obs_cubes = [
            zonal_average.subtract_zonal_average(cube, lat_img, cube_zonal, hdr_zonal)
            for cube, lat_img in zip(obs_cubes, lat_images)
        ]

        data_header = headers[0][1]
        with tools.ignore_warnings(
            'invalid value encountered in log', 'divide by zero encountered'
        ):
            cube_zonal_tb = tools.brightness_temperature(cube_zonal, data_header)
        tb_cubes_maybe_diff = [
            zonal_average.subtract_zonal_average(tb, lat_img, cube_zonal_tb, hdr_zonal)
            for tb, lat_img in zip(tb_cubes, lat_images)
        ]

        obs_cmap = OBS_CMAP_DIFF
        tb_cmap = TB_CMAP
        tb_cmap_maybe_diff = TB_CMAP_DIFF
        difference = True
    else:
        tb_cubes_maybe_diff = tb_cubes
        obs_cmap = OBS_CMAP
        tb_cmap = TB_CMAP
        tb_cmap_maybe_diff = TB_CMAP
        difference = False

    with mstyle.context('dark_background'):  #  type: ignore
        gridspec_kw = dict(
            left=0.04,
            top=0.95,
            right=0.925,
            bottom=0.3,
            width_ratios=[1, 1, 1, 2.3],
            wspace=0.1,
            hspace=0.1,
        )
        sp_rect = [0.075, 0.05, 0.875, 0.205]
        if difference:
            cbar_obs_rect = [0.94, 0.775, 0.01, 0.19]
            cbar_tb_rect = [0.94, 0.5325, 0.01, 0.19]
            cbar_tb_rect_diff = [0.94, 0.29, 0.01, 0.19]
        else:
            cbar_obs_rect = [0.94, 0.65, 0.01, 0.275]
            cbar_tb_rect_diff = []
            cbar_tb_rect = [0.94, 0.325, 0.01, 0.275]

        fig, axs = plt.subplots(
            nrows=len(observations),
            ncols=4,
            squeeze=False,
            gridspec_kw=gridspec_kw,
            facecolor=FIGURE_BACKGROUND_COLOR,
            subplot_kw=dict(facecolor=AXES_BACKGROUND_COLOR),
            figsize=(12, 12),
            dpi=dpi,
        )

        extend = 'both' if clip_percentile else 'neither'
        cbar_obs = fig.colorbar(
            ScalarMappable(cmap=obs_cmap),
            cax=fig.add_axes(cbar_obs_rect),
            extend=extend,
            label='Measured flux{}, MJy/sr'.format(' difference' if difference else ''),
        )
        cbar_tb = fig.colorbar(
            ScalarMappable(cmap=tb_cmap),
            cax=fig.add_axes(cbar_tb_rect),
            extend=extend,
            label='Brightness temperature, K',
        )
        if difference:
            cbar_tb_diff = fig.colorbar(
                ScalarMappable(cmap=tb_cmap_maybe_diff),
                cax=fig.add_axes(cbar_tb_rect_diff),
                extend=extend,
                label='Brightness temperature difference, K',
            )
        else:
            cbar_tb_diff = cbar_tb
        for cbar in [cbar_obs, cbar_tb, cbar_tb_diff]:
            cbar.ax.yaxis.set_label_position('left')
            cbar.ax.yaxis.set_tick_params(rotation=90)

        ax_sp = fig.add_axes(sp_rect)
        ax_sp.plot(wavelengths, spectrum, linewidth=0.5, color='w')
        ax_sp.set_xlabel('Wavelength, µm')
        ax_sp.set_ylabel('Brightness temperature, K')
        ax_sp.set_xlim(min(wavelengths), max(wavelengths))
        for spines in ['top', 'right']:
            ax_sp.spines[spines].set_visible(False)

        fig.suptitle(title, y=0.99)

        # PREPARE SUBPLOTS
        title_kw = dict(size='small')
        map_grid_step = 30
        map_tick_step = 90
        map_grid_kw = dict(
            color=WIREFRAME_COLOR, zorder=0, linestyle=':', linewidth=1, alpha=0.5
        )

        km_xy_lim = []
        for obs in observations:
            nx, ny = obs.get_img_size()
            coords = [
                [-0.5, -0.5],
                [-0.5, ny - 0.5],
                [nx - 0.5, -0.5],
                [nx - 0.5, ny - 0.5],
            ]
            km_xy_lim.extend([obs.radec2km(*obs.xy2radec(x, y)) for x, y in coords])
        km_x_lim, km_y_lim = zip(*km_xy_lim)
        km_x_lim = (min(km_x_lim), max(km_x_lim))
        km_y_lim = (min(km_y_lim), max(km_y_lim))

        for obs_idx, (obs, (primary_header, data_header)) in enumerate(
            zip(observations, headers)
        ):
            ax = axs[obs_idx][0]
            if obs_idx == 0:
                ax.set_title('IFU frame', **title_kw)
            ax.set_ylabel(
                'Dither {}\n{}'.format(
                    primary_header['PATT_NUM'],
                    primary_header['DATE-BEG'][:19].replace('T', ' '),  #  type: ignore
                ),
                **title_kw,
            )
            ax.set_xticks([])
            ax.set_yticks([])

            for ax in (axs[obs_idx][1], axs[obs_idx][2]):
                obs.plot_wireframe_km(
                    ax,
                    aspect_adjustable='datalim',
                    add_title=False,
                    add_axis_labels=False,
                    color=WIREFRAME_COLOR,
                    label_poles=False,
                    zorder=0,
                )
                ax.set_xticks([])
                ax.set_yticks([])

            ax = axs[obs_idx][1]
            if obs_idx == 0:
                ax.set_title(
                    primary_header['TARGNAME'].capitalize()  #  type: ignore
                    + (' aligned' if difference else ' context'),
                    **title_kw,
                )
            if difference:
                ax.set_xlim(*km_x_lim)
                ax.set_ylim(*km_y_lim)
            else:
                ax.set_xlim(-obs.r_eq * radius_factor, obs.r_eq * radius_factor)
                ax.set_ylim(-obs.r_eq * radius_factor, obs.r_eq * radius_factor)

            ax = axs[obs_idx][2]
            if obs_idx == 0:
                ax.set_title(
                    primary_header['TARGNAME'].capitalize()  #  type: ignore
                    + ' aligned',
                    **title_kw,
                )
            ax.set_xlim(*km_x_lim)
            ax.set_ylim(*km_y_lim)

            ax = axs[obs_idx][3]
            if obs_idx == 0:
                ax.set_title('Mapped', **title_kw)
            ax.set_aspect(1, adjustable='box')
            xt = np.arange(0, 360.1, map_tick_step)
            xg = np.arange(0, 360.1, map_grid_step)
            ax.set_xticks(xt)
            if obs_idx == len(observations) - 1:
                ax.set_xticklabels(
                    [
                        '' if x % 90 else f'{x:.0f}°{obs.positive_longitude_direction}'
                        for x in xt
                    ],
                    size='x-small',
                )
            else:
                ax.set_xticklabels(['' for x in xt])
            for x in xg[1:-1]:
                ax.axvline(x, **map_grid_kw)

            yt = np.arange(-90, 90.1, map_tick_step)
            yg = np.arange(-90, 90.1, map_grid_step)
            ax.set_yticks(yt)
            ax.set_yticklabels(
                [
                    ''
                    if y % 90
                    else f'{abs(y):.0f}°{"N" if y > 0 else ("S" if y < 0 else "")}'
                    for y in yt
                ],
                size='x-small',
            )
            for y in yg[1:-1]:
                ax.axhline(y, **map_grid_kw)
            if obs.positive_longitude_direction == 'W':
                ax.set_xlim(360, 0)
            else:
                ax.set_xlim(0, 360)

        handles = []

        def animate_frame(wl_idx: int):
            while handles:
                handles.pop().remove()

            wl_slice = slice(wl_idx, wl_idx + wavelength_step)
            with tools.ignore_warnings('Mean of empty slice'):
                obs_images = [np.nanmean(oc[wl_slice], axis=0) for oc in obs_cubes]
                tb_images = [np.nanmean(tbc[wl_slice], axis=0) for tbc in tb_cubes]
                if difference:
                    tb_images_maybe_diff = [
                        np.nanmean(tbc[wl_slice], axis=0) for tbc in tb_cubes_maybe_diff
                    ]
                else:
                    tb_images_maybe_diff = tb_images

            selected_wavelengths = wavelengths[wl_slice]
            w = float(np.mean(selected_wavelengths))
            s = f'{w:.2f}µm' if wavelength_step > 1 else f'{w:.3f}µm'
            handles.append(
                ax_sp.annotate(
                    s,
                    (w, 1.01),
                    xycoords=('data', 'axes fraction'),  #  type: ignore
                    color='w',
                    ha='center',
                    va='bottom',
                    size='large',
                )
            )
            if len(selected_wavelengths) == 1:
                handles.append(ax_sp.axvline(w, color='w', alpha=0.5, linewidth=0.5))
            else:
                # Add offset of half a wavelength step to prevent it appearing like
                # wavelengths are skipped
                offset = np.mean(np.diff(selected_wavelengths)) / 2
                handles.append(
                    ax_sp.axvspan(
                        min(selected_wavelengths) - offset,
                        max(selected_wavelengths) + offset,
                        color='0.2',
                        zorder=0,
                        linewidth=0,
                    )
                )

            tb_maps = [
                obs.map_img(
                    tbi,
                    degree_interval=map_degree_interval,
                    interpolation=map_interpolation,
                )
                for tbi, obs in zip(tb_images_maybe_diff, observations)
            ]

            lim = get_limits(
                obs_images,
                clip_percentile,
                clip_allowed_range,
                symmetric=difference,
            )
            obs_kwargs = dict(cmap=obs_cmap, vmin=lim[0], vmax=lim[1])
            cbar_obs.set_ticks(
                [0, 1], labels=[format(x, '+.1e' if difference else '.1e') for x in lim]
            )

            lim = get_limits(
                tb_images,
                clip_percentile,
                clip_allowed_range,
                mask_images=mask_images,
                symmetric=False,
            )
            tb_kwargs = dict(cmap=tb_cmap, vmin=lim[0], vmax=lim[1])
            cbar_tb.set_ticks([0, 1], labels=[format(x, '.0f') for x in lim])

            if difference:
                lim = get_limits(
                    tb_images_maybe_diff,
                    clip_percentile,
                    clip_allowed_range,
                    mask_images=mask_images,
                    symmetric=difference,
                )
                tb_kwargs_maybe_diff = dict(
                    cmap=tb_cmap_maybe_diff, vmin=lim[0], vmax=lim[1]
                )
                cbar_tb_diff.set_ticks([0, 1], labels=[format(x, '+.2f') for x in lim])
            else:
                tb_kwargs_maybe_diff = tb_kwargs

            for cbar in [cbar_obs, cbar_tb, cbar_tb_diff]:
                cbar.ax.yaxis.set_label_position('left')
                for idx, l in enumerate(cbar.ax.yaxis.get_ticklabels()):
                    l.set_rotation(90)
                    if clip_percentile:
                        l.set_verticalalignment('center')
                    else:
                        l.set_verticalalignment('top' if idx else 'bottom')

            for obs_idx, (obs, obs_img, tb_img, tb_img_maybe_diff, tb_map) in enumerate(
                zip(observations, obs_images, tb_images, tb_images_maybe_diff, tb_maps)
            ):
                ax = axs[obs_idx][0]
                handles.append(ax.imshow(obs_img, origin='lower', **obs_kwargs))

                ax = axs[obs_idx][1]
                handles.append(
                    ax.pcolormesh(
                        tb_img if difference else obs_img,
                        transform=obs.matplotlib_xy2km_transform(ax),
                        **(tb_kwargs if difference else obs_kwargs),
                    )
                )

                ax = axs[obs_idx][2]
                handles.append(
                    ax.pcolormesh(
                        tb_img_maybe_diff,
                        transform=obs.matplotlib_xy2km_transform(ax),
                        **tb_kwargs_maybe_diff,
                    )
                )

                ax = axs[obs_idx][3]
                if obs.positive_longitude_direction == 'W':
                    extent = [360, 0, -90, 90]
                else:
                    extent = [0, 360, -90, 90]
                handles.append(
                    ax.imshow(
                        tb_map,
                        origin='lower',
                        extent=extent,
                        zorder=1,
                        **tb_kwargs_maybe_diff,
                    )
                )
            # return handles

        animate_frame(0)  # prepare first frame to check everything is working

        tools.check_path(output_path)
        temp_path = output_path + '_animating.mp4'
        idxs = range(0, len(wavelengths), wavelength_step)
        if progress_bar:
            idxs = tqdm.tqdm(idxs, desc='Animating')
        ani = mpanim.FuncAnimation(
            fig,
            animate_frame,
            idxs,
            interval=1000 / fps,
            blit=False,
        )
        writer = mpanim.FFMpegWriter(fps=fps)
        ani.save(
            temp_path,
            savefig_kwargs=dict(facecolor=FIGURE_BACKGROUND_COLOR),
            writer=writer,
            dpi=dpi,
        )
        if progress_bar:
            idxs.update()  # type: ignore
            idxs.close()  # type: ignore # Manually close progress bar

        os.replace(temp_path, output_path)
        if print_output:
            print('Animation {} saved'.format(os.path.split(output_path)[1]))
        plt.close(fig)


def get_limits(
    images: list[np.ndarray],
    clip_percentile: float,
    allowed_range: float,
    mask_images: list[np.ndarray] | None = None,
    symmetric: bool = False,
) -> list[float]:
    if mask_images and np.any(mask_images):
        values = np.hstack([i[m] for i, m in zip(images, mask_images)]).ravel()
    else:
        values = np.array(images).ravel()

    vmin = np.nanmin(values)
    vmax = np.nanmax(values)

    pmin = np.nanpercentile(values, clip_percentile)
    pmax = np.nanpercentile(values, 100 - clip_percentile)
    p_range = pmax - pmin

    if vmax > pmax + allowed_range * p_range:
        vmax = pmax
    if vmin < pmin - allowed_range * p_range:
        vmin = pmin

    if symmetric:
        vm = max(abs(vmin), abs(vmax))  #  type: ignore
        vmin = -vm
        vmax = vm

    return [float(vmin), float(vmax)]


def main():
    parser = argparse.ArgumentParser(
        description='Animate observed FITS cubes',
        epilog=CLI_EXAMPLES,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        argument_default=argparse.SUPPRESS,
    )
    parser.add_argument(
        'input_paths',
        type=str,
        nargs='+',
        help='Paths of JWST FITS cubes to animate',
    )
    parser.add_argument(
        'output_path',
        type=str,
        help='Path of output animation',
    )
    parser.add_argument(
        '--fps',
        type=int,
        help='Frames per second of output animation',
    )
    parser.add_argument(
        '--wavelength_step',
        type=int,
        help='Number of wavelengths averaged per frame of animation',
    )
    parser.add_argument(
        '--map_degree_interval',
        type=float,
        help='Degree interval of plotted map',
    )
    parser.add_argument(
        '--map_interpolation',
        type=str,
        choices=['nearest', 'linear', 'quadratic', 'cubic'],
        help='Interpolation used when projecting mapped data',
    )
    parser.add_argument(
        '--clip_percentile',
        type=float,
        help='Percentile of to set colour limits (to identify outliers)',
    )
    parser.add_argument(
        '--clip_allowed_range',
        type=float,
        help='Tolerance factor for identifying outliers outside the percentile range',
    )
    parser.add_argument(
        '--radius_factor',
        type=float,
        help='Size of context plot limits in multiples of the target\'s equatorial diameter',
    )
    parser.add_argument(
        '--dpi', type=int, help='Resolution (dots per inch) of output animation'
    )
    parser.add_argument(
        '--title_suffix',
        type=str,
        help='Additional suffix to add to the title of the animation',
    )
    parser.add_argument(
        '--zonal_average_path',
        type=str,
        help='Path of zonal average data to plot deviations',
    )
    args = vars(parser.parse_args())
    make_animation(**args)


if __name__ == '__main__':
    main()
