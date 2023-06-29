#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to create useful summary plots of JWST NIRSPEC and MIRI observations.

The summary plots include:
- approximate geometry information
- images at multiple wavelengths
- average spectra
- saturation information

The input files to this script must be navigated using the JWSTSolarSystemPointing code:
https://github.com/JWSTGiantPlanets/NIRSPEC-Toolkit/blob/main/navigate_jwst_observations.py


Requirements
------------
- numpy
- matplotlib
- astopy
- planetmapper


Example Usage
-------------
This script can be run in python, for example: ::

    from jwst_summary_plots import make_summary_plot

    make_summary_plot(
        'path/to/fits/data_s3d.fits',
        'path/to/save/plot.png',
    )

This can also be run directly from the command line by specifying the input FITS and 
output image paths: ::

    python3 jwst_summary_plots path/to/fits/data_s3d.fits path/to/save/plot.png


For more detail, see the docstring of `make_summary_plot()`

See also
--------
`make_miri_saturn_summary_plots.py` and `make_nirspec_summary_plots.py` are examples of
generating plots in bulk.

`jwst_summary_animation.py` creates animations to summarise multiple cubes.

SPICE Kernels
-------------
If you are getting errors from SPICE, it is likely that SPICE kernels aren't being 
loaded correctly when attempting to calculate the coordinate transformations. You
probably want to change the paths defined in PLANETMAPPER_KW below to fix this, or
follow the instructions at 
https://planetmapper.readthedocs.io/en/latest/spice_kernels.html#customising-the-kernel-directory
"""
import os
import pathlib
import sys
import warnings

import matplotlib.pyplot as plt
import matplotlib.scale
import matplotlib.ticker
import numpy as np
import planetmapper
from astropy.io import fits

PLANETMAPPER_KW = {}
if sys.platform == 'linux' and 'PLANETMAPPER_KERNEL_PATH' not in os.environ:
    # Try to set the kernel path automatically when running on ALICE
    PLANETMAPPER_KW['kernel_path'] = '/data/nemesis/jwst/scripts/kernels'

# Customise the plot
DATA_CMAP = 'plasma'
SATURATION_CMAP = 'YlOrRd'
SATURATION_COLOR = 'r'


def main(p_in, p_out):
    make_summary_plot(p_in, p_out)


def make_summary_plot(
    path_in: str,
    path_out: str | None = None,
    show: bool = False,
    vmin_percentile: float = 0,
    vmax_percentile: float = 100,
    plot_brightest_spectrum: bool = True,
):
    """
    Create and save a summary plot of a JWST observation.

    The input FITS file should be navigated and include RA and DEC backplane
    information. The navigation information can be generated using
    https://github.com/JWSTGiantPlanets/NIRSPEC-Toolkit/blob/main/navigate_jwst_observations.py

    If you just want to view the summary plot (rather than saving it), then call
    this function with `show = True`.

    For NIRSPEC observations (which seem to contain spurious hot/dark pixels), this
    function attempts to automatically choose useful limits so that the actual data is
    visible.

    Args:
        path_in: Path of input FITS file.
        path_out: Path to save summary plot image. If `path_out` is `None`, then the
            output path will be automatically generated by replacing `'.fits'` with
            `'.png'` in the input path. If the output directory does not exist, it will
            automatically be created.
        show: Toggle to optionally show the plot instead of saving it.
    """
    fig = plt.figure(figsize=(10, 10), dpi=200)

    ncols = 6
    grid = (6, ncols)
    ax_context = plt.subplot2grid(grid, (0, 0), colspan=ncols, rowspan=2, fig=fig)
    ax_sp = plt.subplot2grid(grid, (2, 0), colspan=ncols, rowspan=2, fig=fig)
    data_row = 4
    sat_row = 5

    with fits.open(path_in) as hdul:
        data = hdul['SCI'].data  # type: ignore
        primary_header = hdul['PRIMARY'].header  # type: ignore
        data_header = hdul['SCI'].header  # type: ignore
        ra_img = hdul['RA'].data  # type: ignore
        dec_img = hdul['DEC'].data  # type: ignore
        basic_navigation = primary_header.get('HIERARCH NAV BASIC', False)
        if basic_navigation:
            lon_img = np.full_like(ra_img, np.nan)
        else:
            lon_img = hdul['LON'].data  # type: ignore
        reduction_notes = get_header_reduction_notes(hdul)
    instrument = primary_header['INSTRUME']

    wl = get_wavelengths(data_header)
    wavelengths = np.linspace(max(wl), min(wl), 20)

    with ignore_warnings('Mean of empty slice'):
        if instrument == 'NIRSPEC':
            sp_fn = np.nanmedian
            sp_type = 'Median'
        else:
            sp_fn = np.nanmean
            sp_type = 'Mean'
    sp_saturated = np.sum(data == 0, axis=(1, 2)) + np.sum(np.isnan(data), axis=(1, 2))

    if instrument == 'NIRSPEC':
        img = np.nanmedian(data, axis=0)
    else:
        img = np.nansum(data, axis=0)
    img[img == 0] = np.nan

    ax = ax_context
    annotate_parts = [
        primary_header['DATE-BEG'],
    ]

    if basic_navigation:
        body = None
    else:
        body = planetmapper.Observation(
            path_in, aberration_correction='CN', **PLANETMAPPER_KW
        )
        body.disc_from_wcs(True, False)

    if body is None:
        planetmapper.utils.format_radec_axes(ax, dec=np.nanmean(dec_img))
    else:
        if body.get_r0() < 1:
            body.plot_wireframe_radec(ax, zorder=-1)
            ax.scatter(
                body.target_ra,
                body.target_dec,
                color='k',
                marker='x',  # type: ignore
                zorder=10,
            )
        else:
            body.plot_wireframe_radec(ax)

        if body.get_r0() > max(img.shape) and not np.all(np.isnan(lon_img)):
            # If the target is bigger than the FOV (and the target is in the FOV), fix the
            # top subplot to be centred on the target
            ax.autoscale_view()
            ax.autoscale(False)
        annotate_parts.append(
            'Sub-obs point: ({ew:.0f}°{ew_s}, {n:.0f}°N)'.format(
                ew=body.subpoint_lon,
                n=body.subpoint_lat,
                ew_s=body.positive_longitude_direction,
            )
        )
        annotate_parts.append(
            'Disc diameter: ~{px:.2f}px ({arcsec:.2f}″)'.format(
                px=body.get_r0() * 2,
                arcsec=body.target_diameter_arcsec,
            )
        )

    ax.annotate(
        '\n'.join(annotate_parts),
        (0.005, 0.99),
        xycoords='axes fraction',
        ha='left',
        va='top',
        size='small',
    )
    ax.annotate(
        '\n'.join(reduction_notes),
        (0.995, 0.99),
        xycoords='axes fraction',
        ha='right',
        va='top',
        size='small',
    )

    ra_offset = primary_header.get('HIERARCH NAV RA_OFFSET', 0)
    dec_offset = primary_header.get('HIERARCH NAV DEC_OFFSET', 0)
    if ra_offset or dec_offset:
        s = f'Navigated using manual offsets:\nRA={ra_offset:+g}″, Dec={dec_offset:+g}″'
    else:
        s = 'Navigated using JWST pointing data\n(typically accurate to ~0.5″)'
    ax.annotate(
        s,
        (0.005, 0.01),
        xycoords='axes fraction',
        ha='left',
        va='bottom',
        size='small',
    )

    vmin_vmax_notes = []
    if vmin_percentile != 0 or vmax_percentile != 100:
        vmin_vmax_notes.append('Image percentile limits')
        vmin_vmax_notes.append(f'min={vmin_percentile}%, max={vmax_percentile}%')
    ax.annotate(
        '\n'.join(vmin_vmax_notes),
        (0.995, 0.01),
        xycoords='axes fraction',
        ha='right',
        va='bottom',
        size='small',
    )

    # ax.contourf(ra_img % 360, dec_img, img, cmap=DATA_CMAP, levels=255)
    ax.pcolormesh(
        ra_img % 360,
        dec_img,
        img,
        cmap=DATA_CMAP,
        vmin=np.nanpercentile(img, vmin_percentile),
        vmax=np.nanpercentile(img, vmax_percentile),
    )

    if 'dither-combined' in primary_header['ASNTABLE']:
        # Dither combinations still have PATT_NUM = 1, so use asn filename as flag
        dither = 'combined'
    else:
        dither = primary_header['PATT_NUM']

    title_parts = [
        primary_header['TARGPROP'],
        f'Dither {dither}',
        primary_header['INSTRUME'],
    ]
    if instrument == 'MIRI':
        title_parts.append(
            'Channel {} {}'.format(primary_header['CHANNEL'], primary_header['BAND'])
        )
        if primary_header['S_RESFRI'] == 'COMPLETE':
            title_parts.append('Residual fringe corrected')
    if instrument == 'NIRSPEC':
        title_parts.append(
            '{} {} {}'.format(
                primary_header['DETECTOR'],
                primary_header['FILTER'],
                primary_header['GRATING'],
            )
        )
    ax.set_title('  |  '.join(title_parts))

    ax = ax_sp
    sp = sp_fn(data, axis=(1, 2))
    data_to_use = img > np.nanpercentile(img, 90)
    sp_max = sp_fn(data[:, data_to_use], axis=1)

    if plot_brightest_spectrum:
        ax.plot(wl, sp_max, color='tab:blue', linewidth=1)
    ax.plot(wl, sp, color='k', linewidth=1)
    if plot_brightest_spectrum:
        for s, c in [
            ('Whole cube', 'k'),
            (' ' * 20 + 'Brightest 10% spaxels', 'tab:blue'),
        ]:
            ax.annotate(
                s,
                (0.005, 0.99),
                xycoords='axes fraction',
                ha='left',
                va='top',
                size='x-small',
                color=c,
                zorder=0,
            )

    ax.set_xlim(min(wl), max(wl))
    ax.set_xlabel('Wavelength (µm)')

    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    sp_units = data_header['BUNIT']
    ax.set_ylabel(f'{sp_type} spectrum ({sp_units})')
    if np.sum(sp_max < 0) > len(sp_max) / 3 or np.all(~np.isfinite(sp_max)):
        ax.set_ylabel(ax.get_ylabel() + ' (linear scale)')
    else:
        ax.set_yscale('log')
    if instrument == 'NIRSPEC':
        ylim = ax.get_ylim()
        ymin = max(
            ylim[0], min(max(np.nanpercentile(sp, 5) / 10, 5e-2), np.nanmin(sp_max))
        )
        if ymin > 0 or ax.get_yscale() != 'log':
            ax.set_ylim(bottom=ymin)

    ax = ax.twinx()
    ax.plot(wl, sp_saturated, color=SATURATION_COLOR, linewidth=0.5)
    ax.set_ylim(*ax.get_ylim())
    ax.fill_between(wl, sp_saturated, color=SATURATION_COLOR, linewidth=0, alpha=0.1)
    ax.set_ylim(top=data.shape[1] * data.shape[2])
    ax.set_ylabel('Invalid or fully saturated pixels')
    ax.spines['right'].set_color(SATURATION_COLOR)
    ax.yaxis.label.set_color(SATURATION_COLOR)
    ax.tick_params(axis='y', colors=SATURATION_COLOR)

    wavelengths = np.linspace(min(wl), max(wl), ncols + 1)
    for wl_idx in range(ncols):
        ax_data = plt.subplot2grid(grid, (data_row, wl_idx), fig=fig)
        ax_sat = plt.subplot2grid(grid, (sat_row, wl_idx), fig=fig)
        wl0 = wavelengths[wl_idx]
        wl1 = wavelengths[wl_idx + 1]
        img_data = data[np.logical_and(wl >= wl0, wl < wl1)]
        if instrument == 'NIRSPEC':
            img_data = np.nanmedian(img_data, axis=0)
        else:
            img_data = np.nansum(img_data, axis=0)
        img_sat = np.nanmean(
            data[np.logical_and(wl >= wl0, wl < wl1)] == 0, axis=0
        ) + np.nanmean(np.isnan(data[np.logical_and(wl >= wl0, wl < wl1)]), axis=0)
        img_sat[img_sat == 0] = np.nan

        ax_data.imshow(
            img_data,
            cmap=DATA_CMAP,
            origin='lower',
            vmin=np.nanpercentile(img_data, vmin_percentile),
            vmax=np.nanpercentile(img_data, vmax_percentile),
        )
        ax_sat.imshow(img_sat, cmap=SATURATION_CMAP, origin='lower', vmin=-0.15, vmax=1)

        for ax in [ax_data, ax_sat]:
            ax.set_xticks([])
            ax.set_yticks([])

        ax_data.set_xlabel(f'{wl0:.2f} - {wl1:.2f}µm')

        if wl_idx % 2:
            ax_sp.axvspan(wl0, wl1, color='k', zorder=0, alpha=0.05, linewidth=0)

        if wl_idx == 0:
            ax = ax_sat
            if body is not None:
                if body.get_r0() > 2:
                    body.plot_wireframe_xy(ax)
                else:
                    ax.scatter(
                        body.get_x0(),
                        body.get_y0(),
                        color='k',
                        marker='x',  # type: ignore
                        zorder=5,
                    )

                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title('')
                rot = body.get_rotation() + body.north_pole_angle()

                length = 0.75
                dx = np.sin(np.deg2rad(rot)) * length
                dy = np.cos(np.deg2rad(rot)) * length

                ax.annotate(
                    'N',
                    weight='bold',
                    size='small',
                    xy=(0.5 - dx / 2, 0.5 - dy / 2),
                    xytext=(0.5 + dx / 2, 0.5 + dy / 2),
                    xycoords='axes fraction',
                    arrowprops=dict(arrowstyle='<-', color='b'),
                    color='b',
                )

            ax_data.set_ylabel(
                f'IFU plane\n{"median" if instrument == "NIRSPEC" else "summed"} data'
            )
            ax_sat.set_ylabel('% Invalid\nor saturated')

    fig.tight_layout()
    if show:
        plt.show()
        return

    if path_out is None:
        path_out = path_in.replace('.fits', '.png')
    check_path(path_out)
    fig.savefig(path_out)
    plt.close(fig)


def get_wavelengths(
    header: fits.Header,
    warn=True,
    axis=3,
) -> np.ndarray:
    """
    Get wavelengths from FITS header.
    """
    if header[f'CTYPE{axis}'] != 'WAVE' and warn:
        print(f'WARNING: header item CTYPE{axis} = {header[f"CTYPE{axis}"]}')
    crval3 = header[f'CRVAL{axis}']
    try:
        cd3_3 = header[f'CD{axis}_{axis}']
    except KeyError:
        cd3_3 = header[f'CDELT{axis}']
    naxis3 = header[f'NAXIS{axis}']
    wavl = crval3 + cd3_3 * np.arange(0, naxis3)
    return wavl


def get_header_reduction_notes(
    hdul: fits.HDUList,
    key_prefix: str = 'REDUCT',
) -> list[str]:
    notes = []
    n = 1
    header = hdul['PRIMARY'].header  #  type: ignore
    while True:
        key = f'{key_prefix}{n}'
        if key not in header:
            break
        notes.append(header[key])
        n += 1
    return notes


def check_path(path) -> None:
    """
    Checks if file path's directory tree exists, and creates it if necessary.

    Assumes path is to a file if `os.path.split(path)[1]` contains '.',
    otherwise assumes path is to a directory.

    Parameters
    ----------
    path : str
        Path to directory to check.
    """
    if os.path.isdir(path):
        return
    if '.' in os.path.split(path)[1]:
        path = os.path.split(path)[0]
        if os.path.isdir(path):
            return
    if path == '':
        return
    print('Creating directory path "{}"'.format(path))
    pathlib.Path(path).mkdir(parents=True, exist_ok=True)


class ignore_warnings(warnings.catch_warnings):
    """
    Context manager to hide FITS `Card is too long, comment will be truncated` warnings.
    """

    def __init__(self, *warining_strings: str, **kwargs):
        super().__init__(**kwargs)
        self.warning_strings = warining_strings

    def __enter__(self):
        out = super().__enter__()
        for ws in self.warning_strings:
            warnings.filterwarnings('ignore', ws)
        return out


if __name__ == '__main__':
    main(*sys.argv[1:])  # pylint: disable=too-many-function-args
