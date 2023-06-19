import numpy as np
from astropy import units as u


def save_spx(
    path: str,
    wavelengths: np.ndarray,
    spectrum: np.ndarray,
    error: np.ndarray,
    lon: float,
    lat: float,
    phase: float,
    emission: float,
    azimuth: float,
    pixelarea_arcsecsq: float,
    pixelarea_steradians: float,
    fwhm: float = 0,
) -> None:
    # Based on Henrik's solar system pointing code
    spectrum = convert_MJysr_to_Wcm2srum(
        wavelengths, spectrum, pixelarea_arcsecsq, pixelarea_steradians
    )
    error = convert_MJysr_to_Wcm2srum(
        wavelengths, error, pixelarea_arcsecsq, pixelarea_steradians
    )

    ngeom = 1
    nconv = wavelengths.shape[0]
    nav = 1
    wgeom = 1

    # Generate the header
    header = []
    header.append([fwhm, lat, lon, ngeom])
    header.append([nconv])
    header.append([nav])
    header.append([lat, lon, phase, emission, azimuth, wgeom])

    with open(path, 'w', encoding='utf-8') as f:
        for line in header:
            f.write('\t'.join('{:.3f}'.format(x) for x in line))
            f.write('\n')
        for i, w in np.ndenumerate(wavelengths):
            line = [wavelengths[i], spectrum[i], error[i]]
            f.write('\t'.join('{:.6e}'.format(x) for x in line))
            f.write('\n')


def convert_MJysr_to_Wcm2srum(
    wave, spx_MJysr, pixelarea_arcsecsq, pixelarea_steradians
):
    # Unit converstion from MJy/sr to W/cm2/sr/micron as per Leigh Fletcher
    # From Henrik's solar system pointing code

    # Add the correct units
    spx_MJysr = spx_MJysr * u.MJy / u.sr

    # Convert "per sr" to "per square arcsecond"
    spx_Jy_per_arcsec = spx_MJysr.to(u.Jy / (u.arcsec * u.arcsec))

    # Integrate over solid angle to get irradiance (Jy)
    spx_Jy = spx_Jy_per_arcsec * pixelarea_arcsecsq * u.arcsec * u.arcsec
    # spx_Jy=spx_MJysr * 1e-6 * PIXAR_SR*u.sr

    # Now convert Jy to  W/cm2/sr/cm-1
    c = 2.99792458e10  # *u.cm/u.s  #cm/s
    corrn = 1e-26  # *u.Jy/(u.W/(u.m*u.m)/u.Hz) #(Jy -> W/m2/Hz)
    corrn = corrn / (1.0e4)  # W/m2/Hz -> W/cm2/Hz
    corrn = corrn * c  # W/cm2/Hz -> W/cm2/cm-1
    corrn = corrn / pixelarea_steradians  # *u.sr    # W/cm2/cm-1 -> W/cm2/sr/cm-1

    spx_Wcm2srcm = spx_Jy * corrn

    # Convert  W/cm2/sr/cm-1 to TB

    h = 6.626e-34
    c = 2.9979e8
    k = 1.3806e-23
    spec = spx_Wcm2srcm * 1e4 / 100.0 / u.Jy  # W/cm2/sr/cm-1 -> W/m2/sr/m-1
    v = 100.0 * 1e4 / wave  # wavenumber in m-1
    c1 = 2 * h * c * c
    c2 = h * c / k
    a = c1 * v * v * v / spec
    TB = c2 * v / (np.log(a + 1))

    # Now convert back to W/cm2/sr/µm
    l = wave * 1e-6
    a = 2 * h * c * c / (l**5)
    b = np.exp(h * c / (l * k * TB)) - 1
    spx_Wcm2srum = (a / b) / 1e4 / 1e6

    # return spx_Jy_per_arcsec, spx_Jy, spx_Wcm2srcm, spx_Wcm2srum, TB
    return spx_Wcm2srum
