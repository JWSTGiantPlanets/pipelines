"""General useful tools"""
import datetime
import hashlib
import itertools
import math
import os
import pathlib
import pickle
import re
import subprocess
import warnings
from typing import Any, Callable, ParamSpec, TypeVar

import astropy.units as u
import numpy as np
from astropy.io import fits

import function_cache

P = ParamSpec('P')
T = TypeVar('T')


# Warnings
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


# Maths
def nearest_idx(array, value: float) -> int:
    """
    Find nearest index in array to a value.

    More general than list.index() method which only works for an exact match.

    >>> tools.maths.nearest_idx([1,2,3], 2.1)
    2

    Parameters
    ----------
    array : list

    value : float or int
        Value to find closest index for

    Returns
    -------
    int
    """
    diff_arr = np.asarray(np.abs(np.array(array) - value), dtype=float)
    diff_arr[np.isnan(diff_arr)] = np.inf
    return diff_arr.argmin()  #  type: ignore


def nearest_idx_nd(array, values: np.ndarray) -> np.ndarray:
    out = np.full(values.shape, -1, dtype=int)
    for index, value in np.ndenumerate(values):
        if np.isfinite(value):
            out[index] = nearest_idx(array, value)
    return out


def clip_percentile(data, bottom, top=None):
    if top is None:
        top = 100 - bottom
    return np.clip(data, np.nanpercentile(data, bottom), np.nanpercentile(data, top))


def normalise(
    values, percentile: bool | float = False, top=1, bottom=0, single_value=None
):
    """
    Normalise iterable.

    Parameters
    ----------
    values : list
        Iterable to normalise.

    percentile : bool or float or int
        Percentile to set top and bottom values e.g. percentile = 5 sets the
        5th percentile of values as 0 and the 95th as 1.

    top, bottom : float
        Limits of normalised values.

    single_value : float
        Value to return if only one value is present in values.

    Returns
    -------
    Normalised values
    """
    assert top > bottom
    values = np.array(values)
    if single_value is not None and len(set(values)) == 1:
        return np.full(values.shape, single_value)
    if percentile:
        vmin = np.nanpercentile(values, percentile)
        vmax = np.nanpercentile(values, 100 - percentile)
    else:
        vmin = np.nanmin(values)
        vmax = np.nanmax(values)

    # Put into 0 to 1 range
    if vmax != vmin:
        values = (values - vmin) / (vmax - vmin)
    else:
        values = values - vmin
    return values * (top - bottom) + bottom


def rms(a: np.ndarray, b: np.ndarray) -> float:
    c = a - b
    # Use explicit multiplication instead of power as this is slightly faster
    return math.sqrt(np.mean(c * c))


# File IO
def check_path(path):
    """
    Checks if file path's directory tree exists, and creates it if necessary.

    Assumes path is to a file if `os.path.split(path)[1]` contains '.',
    otherwise assumes path is to a directory.

    Parameters
    ----------
    path : str
        Path to directory to check.
    """
    path = os.path.expandvars(os.path.expanduser(path))
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


def get_wavelengths(
    header,
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
    crpix3 = header.get(f'CRPIX{axis}', 1)
    wavl = crval3 + cd3_3 * (np.arange(0, naxis3) + crpix3 - 1)
    return wavl


# Function result caching
def cached(fn: Callable[P, T], *args: P.args, **kwargs: P.kwargs) -> T:
    key = (
        fn.__module__,
        fn.__name__,
        generic_hash(args),
        generic_hash(kwargs),
    )
    if key not in function_cache.CACHE:
        function_cache.CACHE[key] = fn(*args, **kwargs)
    return function_cache.CACHE[key]


def clear_cached_result(
    fn: Callable[P, Any], *args: P.args, **kwargs: P.kwargs
) -> None:
    key = (
        fn.__module__,
        fn.__name__,
        generic_hash(args),
        generic_hash(kwargs),
    )
    function_cache.CACHE.pop(key, None)


def clear_cached_fn(fn: Callable) -> None:
    num = 0
    for k in list(function_cache.CACHE.keys()):
        if k[0] == fn.__module__ and k[1] == fn.__name__:
            function_cache.CACHE.pop(k)
            num += 1
    print(f'Removed {num} cached results')


def clear_cached_all() -> None:
    function_cache.CACHE.clear()


clear_cached = clear_cached_all


def generic_hash(obj):
    h = hashlib.sha1()
    try:
        h.update(pickle.dumps(obj))
    except pickle.PicklingError:
        h.update(repr(obj).encode('utf-8'))
    return h.hexdigest()


def brightness_temperature(
    sp: np.ndarray,
    data_header: fits.Header | None = None,
    *,
    pixar_a2: float | None = None,
    pixar_sr: float | None = None,
    wl: list[float] | np.ndarray | None = None,
) -> np.ndarray:
    if pixar_a2 is None:
        pixar_a2 = data_header['PIXAR_A2']
    if pixar_sr is None:
        pixar_sr = data_header['PIXAR_SR']
    if wl is None:
        wl = get_wavelengths(data_header)

    # Unit converstion from MJy/sr to W/cm2/sr/micron as per Leigh Fletcher

    spx_MJysr = sp * u.MJy / u.sr

    # Convert "per sr" to "per square arcsecond"
    spx_Jy_per_arcsec = spx_MJysr.to(u.Jy / (u.arcsec * u.arcsec))

    # Integrate over solid angle to get irradiance (Jy)
    spx_Jy = spx_Jy_per_arcsec * pixar_a2 * u.arcsec * u.arcsec
    # spx_Jy=spx_MJysr * 1e-6 * PIXAR_SR*u.sr

    # Now convert Jy to  W/cm2/sr/cm-1
    c = 2.99792458e10  # *u.cm/u.s  #cm/s
    corrn = 1e-26  # *u.Jy/(u.W/(u.m*u.m)/u.Hz) #(Jy -> W/m2/Hz)
    corrn = corrn / (1.0e4)  # W/m2/Hz -> W/cm2/Hz
    corrn = corrn * c  # W/cm2/Hz -> W/cm2/cm-1
    corrn = corrn / pixar_sr  # *u.sr    # W/cm2/cm-1 -> W/cm2/sr/cm-1

    spx_Wcm2srcm = spx_Jy * corrn

    # Convert  W/cm2/sr/cm-1 to TB

    h = 6.626e-34
    c = 2.9979e8
    k = 1.3806e-23
    spec = spx_Wcm2srcm * 1e4 / 100.0 / u.Jy  # W/cm2/sr/cm-1 -> W/m2/sr/m-1
    v = 100.0 * 1e4 / wl  # wavenumber in m-1
    while len(v.shape) < len(sp.shape):
        # work for cubes passed to sp
        v = np.expand_dims(v, -1)
    c1 = 2 * h * c * c
    c2 = h * c / k
    a = c1 * v * v * v / spec
    tb = c2 * v / (np.log(a + 1))
    return np.array(tb, dtype=float)


def add_header_reduction_note(
    hdul: fits.HDUList,
    note: str,
    key_prefix: str = 'REDUCT',
    do_date: bool = True,
) -> None:
    n = 1
    header = hdul['PRIMARY'].header  #  type: ignore
    while True:
        key = f'{key_prefix}{n}'
        if key in header:
            n += 1
            continue
        if n == 1:
            # First reduction note, so add separator
            sep = 'Custom pipeline reductions'
            sep = ' ' * 72 + sep + ' ' * (72 * 2 - len(sep))
            header.append(('', sep), useblanks=False, bottom=True)
        header.append((key, note), useblanks=n > 1, bottom=n == 1)
        if do_date:
            date = datetime.datetime.now().isoformat()
            header.append((f'HIERARCH {key} DATE', date, 'Date reduction perfomed'))
        break


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


class KeepMissingDict(dict):
    def __missing__(self, key) -> str:
        return '{' + key + '}'


def all_combinations(**kwargs) -> list[dict]:
    keys, vals = zip(*[(k, v) for k, v in kwargs.items() if v])
    return [dict(zip(keys, instance)) for instance in itertools.product(*vals)]


def cprint(*msg, fg=None, bg=None, style=None, skip_print=False, sep=' ', **kwargs):
    """
    Prints coloured and formatted text.
    Parameters
    ----------
    msg
        Message to print.
    fg, bg : {'k', 'r', 'g', 'b', 'y', 'm', 'c', 'w'}
        Foreground and background colours (see code for options).
    style : {'b', 'f', 'i', 'u', 'x', 'y', 'r', 'h', 's'}
        Formatting style to apply to text. Can be multiple values, e.g. 'bi'
        for bold and italic style.
    """
    colcode = {
        'k': 0,  # black
        'r': 1,  # red
        'g': 2,  # green
        'y': 3,  # yellow
        'b': 4,  # blue
        'm': 5,  # magenta
        'c': 6,  # cyan
        'w': 7,  # white
    }

    fmtcode = {
        'b': 1,  # bold
        'f': 2,  # faint
        'i': 3,  # italic
        'u': 4,  # underline
        'x': 5,  # blinking
        'y': 6,  # fast blinking
        'r': 7,  # reverse
        'h': 8,  # hide
        's': 9,  # strikethrough
    }

    # Properties
    props = []
    if isinstance(style, str):
        props = [fmtcode[s] for s in style]
    if isinstance(fg, str):
        props.append(30 + colcode[fg])
    if isinstance(bg, str):
        props.append(40 + colcode[bg])

    # Display
    msg = sep.join(
        str(x) for x in msg
    )  # Reproduce print() behaviour for easy translation
    props = ';'.join(str(x) for x in props)

    if props:
        msg = '\x1b[%sm%s\x1b[0m' % (props, msg)

    if not skip_print:
        print(msg, **kwargs)

    return msg


def print_bar_chart(
    labels,
    bars=None,
    formats=None,
    print_values=True,
    max_label_length=None,
    sort=False,
    **kwargs,
):
    """
    Print bar chart of data
    Parameters
    ----------
    labels : array
        Labels of bars, or bar values if `bars is None`.
    bars : array
        List of bar lengths.
    formats : array
        List of bar formats to be passed to `cprint()`.
    print_values : bool
        Toggle printing bar lengths.
    max_label_length : int
        Set length to trim labels, None for no trimming.
    sort : bool
        Toggle sorting of bars by size.
    **kwargs
        Arguments passed to `cprint()` for every bar.
    """
    if bars is None:
        bars = labels.copy()
        labels = ['' for _ in bars]
    bars = list(bars)
    labels = [str(l) for l in labels]
    if max_label_length is None:
        max_label_length = max([len(l) for l in labels] + [1])
    else:
        labels = clip_string_list(labels, max_label_length)
    labels = [f'{l:{max_label_length}s}' for l in labels]
    if sort:
        if formats is not None and not isinstance(formats, str):
            bars, labels, formats = zip(*sorted(zip(bars, labels, formats)))
        else:
            bars, labels = zip(*sorted(zip(bars, labels)))
    if print_values:
        fmt = '.2e'
        if isinstance(print_values, str):
            fmt = print_values
        value_strs = [f'{v:{fmt}}' for v in bars]
        labels = [f'{l}|{v}' for l, v in zip(labels, value_strs)]
    max_label_length = max(len(l) for l in labels)
    max_length = get_console_width() - max_label_length - 2
    for idx, label in enumerate(labels):
        kw = {**kwargs}
        if formats:
            if formats == 'auto':
                if bars[idx] / sum(bars) > 0.5:
                    kw.update(fg='y', style='b')
                elif bars[idx] / sum(bars) > 0.1:
                    kw.update(fg='g', style='b')
                elif bars[idx] / sum(bars) > 0.01:
                    kw.update(fg='b', style='b')
                else:
                    kw.update(fg='w', style='f')
            elif formats == 'extreme':
                if bars[idx] == max(bars):
                    kw.update(fg='g', style='b')
                elif bars[idx] == min(bars):
                    kw.update(fg='r', style='b')
                else:
                    kw.update(fg='b', style='b')
            else:
                kw.update(formats[idx])

        chrs = ' ▏▎▍▌▋▊▉█'
        length = max_length * bars[idx] / max(bars)
        decimal_idx = (length - int(length)) * (len(chrs) - 1)
        decimal_idx = int(np.round(decimal_idx))
        bar = chrs[-1] * int(length) + chrs[decimal_idx]
        bar = bar.rstrip(' ')
        cprint(f'{label:{max_label_length}s}|{bar}', **kw)


def get_console_width(fallback=75, maximum=98):
    """
    Attempts to find console width, otherwise uses fallback provided.
    Parameters
    ----------
    fallback : int
        Default width value if `stty size` fails.
    Returns
    -------
    width : int
        Console width.
    """
    if test_if_ipython():
        return fallback
    try:
        _, width = subprocess.check_output(
            ['stty', 'size'], stderr=subprocess.PIPE
        ).split()
    # pylint: disable-next=bare-except
    except:
        width = fallback
    width = int(width)
    if maximum and width > maximum:
        width = maximum
    return width


def test_if_ipython():
    """Detect if script is running in IPython console"""
    try:
        return __IPYTHON__  # type: ignore
    except NameError:
        return False


def clip_string(s, max_len, continue_str='…'):
    """
    Takes string and clips to certain length if needed.
    Parameters
    ----------
    s : str
        String to clip
    max_len : int
        Maximum allowed string length
    continue_str : str
        String used to indicate string has been clipped
    Returns
    -------
    clipped string
    """
    return s if len(s) <= max_len else s[: max_len - len(continue_str)] + continue_str


def clip_string_list(a, max_len, **kwargs):
    """
    Takes a list of strings and clips them to a certain length if needed.
    Parameters
    ----------
    a : list of str
    Other parameters passed to clip_string()
    Returns
    -------
    clipped list
    """
    return [clip_string(s, max_len, **kwargs) for s in a]


def log(*messages: Any, time: bool = True) -> None:
    """
    Print a message with a timestamp.

    Args:
        *messages: Messages passed to `print()`.
        time: Toggle showing the timestamp.
    """
    prefix = datetime.datetime.now().strftime('%H:%M:%S') if time else ' ' * 8
    print(prefix, *messages, flush=True)


def latex_formula_string(string_in):
    if '\n' in string_in:
        return '\n'.join(latex_formula_string(p) for p in string_in.split('\n'))
    split_list = re.compile(r'(\d+|\.)').split(string_in)
    split_list = [s for s in split_list if s]
    out = ''
    for idx, s in enumerate(split_list):
        if s.isdigit():
            if idx and split_list[idx - 1] != '.':
                s = '_' + s
        if s == '.':
            s = r' \cdot '
            if idx and split_list[idx - 1].isdigit():
                if idx == 1 or split_list[idx - 2] == '.':
                    s = '.'
        out += s
    return r'$\mathrm{' + out + '}$'


def make_rgb_images(
    *images_data: list | np.ndarray, percentile: float = 0
) -> list[np.ndarray]:
    # Get min and max values for each RGB channel
    rgb_values = [
        np.concatenate([img[idx].ravel() for img in images_data]) for idx in range(3)
    ]
    minima = np.nanpercentile(rgb_values, percentile, axis=1)
    maxima = np.nanpercentile(rgb_values, 100 - percentile, axis=1)

    # Change shape to (3, 1, 1) for broadcasting
    minima = minima[:, None, None]
    maxima = maxima[:, None, None]

    # Normalize and convert to RGB
    out = []
    for image_data in images_data:
        image_data = (image_data - minima) / (maxima - minima)
        image_data = np.clip(image_data, 0, 1)
        out.append(np.moveaxis(image_data, 0, -1))
    return out


def merge_nested_dicts(*dicts: dict) -> dict:
    """
    Merge nested dictionaries.

    Args:
        dicts: Dictionaries to merge.
    """
    result = {}
    for d in dicts:
        for key, value in d.items():
            if isinstance(value, dict):
                result[key] = merge_nested_dicts(result.get(key, {}), value)
            else:
                result[key] = value
    return result
