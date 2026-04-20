"""Normalize matched catalog rows into long-form photometry."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from astropy.table import Table, vstack


AB_ZEROPOINT_JY = 3631.0

WISE_VEGA_ZEROPOINT_JY: dict[str, float] = {
    "W1": 309.540,
    "W2": 171.787,
    "W3": 31.674,
    "W4": 8.363,
}

TWOMASS_VEGA_ZEROPOINT_JY: dict[str, float] = {
    "J": 1594.0,
    "H": 1024.0,
    "Ks": 666.7,
}

SPECLITE_NAMES: dict[str, str] = {
    "FUV_galex": "galex-fuv",
    "NUV_galex": "galex-nuv",
    "u_sdss": "sdss2010-u",
    "g_sdss": "sdss2010-g",
    "r_sdss": "sdss2010-r",
    "i_sdss": "sdss2010-i",
    "z_sdss": "sdss2010-z",
    "g_legacy": "decam2014-g",
    "r_legacy": "decam2014-r",
    "z_legacy": "decam2014-z",
    "W1": "wise2010-W1",
    "W2": "wise2010-W2",
    "W3": "wise2010-W3",
    "W4": "wise2010-W4",
    "W1_legacy": "wise2010-W1",
    "W2_legacy": "wise2010-W2",
    "J_2mass": "twomass-J",
    "H_2mass": "twomass-H",
    "Ks_2mass": "twomass-Ks",
    "S9W_akari": "akari-S9W",
    "L18W_akari": "akari-L18W",
    "N60_akari": "akari-N60",
    "WIDE-S_akari": "akari-WIDE-S",
    "WIDE-L_akari": "akari-WIDE-L",
    "N160_akari": "akari-N160",
    "F12_iras": "iras-12",
    "F25_iras": "iras-25",
    "F60_iras": "iras-60",
    "F100_iras": "iras-100",
}

PSF_FWHM_ARCSEC: dict[str, float] = {
    "FUV_galex": 4.3,
    "NUV_galex": 5.3,
    "u_sdss": 1.53,
    "g_sdss": 1.44,
    "r_sdss": 1.32,
    "i_sdss": 1.26,
    "z_sdss": 1.29,
    "g_legacy": 1.30,
    "r_legacy": 1.20,
    "z_legacy": 1.10,
    "W1": 6.08,
    "W2": 6.84,
    "W3": 7.36,
    "W4": 11.99,
    "W1_legacy": 6.08,
    "W2_legacy": 6.84,
    "J_2mass": 2.5,
    "H_2mass": 2.5,
    "Ks_2mass": 2.5,
    "S9W_akari": 9.4,
    "L18W_akari": 10.4,
    "N60_akari": 26.8,
    "WIDE-S_akari": 26.8,
    "WIDE-L_akari": 44.2,
    "N160_akari": 44.2,
    "F12_iras": 30.0,
    "F25_iras": 30.0,
    "F60_iras": 60.0,
    "F100_iras": 120.0,
}


@dataclass(frozen=True)
class BandSpec:
    """Column mapping for one catalog photometric band."""

    band: str
    filter_name: str
    mag_col: str | Sequence[str]
    err_col: str | Sequence[str]
    system: str
    zero_jy: float


@dataclass(frozen=True)
class FluxBandSpec:
    """Column mapping for one catalog flux-density band."""

    band: str
    filter_name: str
    flux_col: str
    err_col: str
    quality_col: str | None = None
    err_is_percent: bool = False


@dataclass(frozen=True)
class NanomaggyBandSpec:
    """Column mapping for one catalog band stored as AB nanomaggies."""

    band: str
    filter_name: str
    flux_col: str | Sequence[str]
    ivar_col: str | Sequence[str]
    psf_col: str | Sequence[str] | None = None


CATALOG_BAND_SPECS: dict[str, tuple[BandSpec, ...]] = {
    "galex_ais": (
        BandSpec("FUV", "FUV_galex", "FUVmag", "e_FUVmag", "ab", AB_ZEROPOINT_JY),
        BandSpec("NUV", "NUV_galex", "NUVmag", "e_NUVmag", "ab", AB_ZEROPOINT_JY),
    ),
    "sdss_dr16": (
        BandSpec("u", "u_sdss", "umag", "e_umag", "ab", AB_ZEROPOINT_JY),
        BandSpec("g", "g_sdss", "gmag", "e_gmag", "ab", AB_ZEROPOINT_JY),
        BandSpec("r", "r_sdss", "rmag", "e_rmag", "ab", AB_ZEROPOINT_JY),
        BandSpec("i", "i_sdss", "imag", "e_imag", "ab", AB_ZEROPOINT_JY),
        BandSpec("z", "z_sdss", "zmag", "e_zmag", "ab", AB_ZEROPOINT_JY),
    ),
    "allwise": (
        BandSpec("W1", "W1", "W1mag", "e_W1mag", "vega", WISE_VEGA_ZEROPOINT_JY["W1"]),
        BandSpec("W2", "W2", "W2mag", "e_W2mag", "vega", WISE_VEGA_ZEROPOINT_JY["W2"]),
        BandSpec("W3", "W3", "W3mag", "e_W3mag", "vega", WISE_VEGA_ZEROPOINT_JY["W3"]),
        BandSpec("W4", "W4", "W4mag", "e_W4mag", "vega", WISE_VEGA_ZEROPOINT_JY["W4"]),
    ),
    "2mass": (
        BandSpec("J", "J_2mass", "Jmag", "e_Jmag", "vega", TWOMASS_VEGA_ZEROPOINT_JY["J"]),
        BandSpec("H", "H_2mass", "Hmag", "e_Hmag", "vega", TWOMASS_VEGA_ZEROPOINT_JY["H"]),
        BandSpec("Ks", "Ks_2mass", "Kmag", "e_Kmag", "vega", TWOMASS_VEGA_ZEROPOINT_JY["Ks"]),
    ),
    "2mass_xsc": (
        BandSpec("J.ext", "J_2mass", ("J.ext", "Jmag"), ("e_J.ext", "e_Jmag"), "vega", TWOMASS_VEGA_ZEROPOINT_JY["J"]),
        BandSpec("H.ext", "H_2mass", ("H.ext", "Hmag"), ("e_H.ext", "e_Hmag"), "vega", TWOMASS_VEGA_ZEROPOINT_JY["H"]),
        BandSpec("K.ext", "Ks_2mass", ("K.ext", "Kmag"), ("e_K.ext", "e_Kmag"), "vega", TWOMASS_VEGA_ZEROPOINT_JY["Ks"]),
    ),
}

CATALOG_FLUX_SPECS: dict[str, tuple[FluxBandSpec, ...]] = {
    "akari_irc": (
        FluxBandSpec("S9W", "S9W_akari", "S09", "e_S09", quality_col="q_S09"),
        FluxBandSpec("L18W", "L18W_akari", "S18", "e_S18", quality_col="q_S18"),
    ),
    "akari_fis": (
        FluxBandSpec("N60", "N60_akari", "S65", "e_S65", quality_col="q_S65"),
        FluxBandSpec("WIDE-S", "WIDE-S_akari", "S90", "e_S90", quality_col="q_S90"),
        FluxBandSpec("WIDE-L", "WIDE-L_akari", "S140", "e_S140", quality_col="q_S140"),
        FluxBandSpec("N160", "N160_akari", "S160", "e_S160", quality_col="q_S160"),
    ),
    "iras_psc": (
        FluxBandSpec("F12", "F12_iras", "Fnu_12", "e_Fnu_12", quality_col="q_Fnu_12", err_is_percent=True),
        FluxBandSpec("F25", "F25_iras", "Fnu_25", "e_Fnu_25", quality_col="q_Fnu_25", err_is_percent=True),
        FluxBandSpec("F60", "F60_iras", "Fnu_60", "e_Fnu_60", quality_col="q_Fnu_60", err_is_percent=True),
        FluxBandSpec("F100", "F100_iras", "Fnu_100", "e_Fnu_100", quality_col="q_Fnu_100", err_is_percent=True),
    ),
}

_LEGACY_NANOMAGGY_SPECS = (
    NanomaggyBandSpec(
        "g",
        "g_legacy",
        ("FLUX_G", "flux_g"),
        ("FLUX_IVAR_G", "flux_ivar_g"),
        ("PSFSIZE_G", "psfsize_g"),
    ),
    NanomaggyBandSpec(
        "r",
        "r_legacy",
        ("FLUX_R", "flux_r"),
        ("FLUX_IVAR_R", "flux_ivar_r"),
        ("PSFSIZE_R", "psfsize_r"),
    ),
    NanomaggyBandSpec(
        "z",
        "z_legacy",
        ("FLUX_Z", "flux_z"),
        ("FLUX_IVAR_Z", "flux_ivar_z"),
        ("PSFSIZE_Z", "psfsize_z"),
    ),
    NanomaggyBandSpec(
        "W1",
        "W1_legacy",
        ("FLUX_W1", "flux_w1"),
        ("FLUX_IVAR_W1", "flux_ivar_w1"),
    ),
    NanomaggyBandSpec(
        "W2",
        "W2_legacy",
        ("FLUX_W2", "flux_w2"),
        ("FLUX_IVAR_W2", "flux_ivar_w2"),
    ),
)

CATALOG_NANOMAGGY_SPECS: dict[str, tuple[NanomaggyBandSpec, ...]] = {
    "desi_legacy_dr8_north": _LEGACY_NANOMAGGY_SPECS,
    "desi_legacy_dr8_south": _LEGACY_NANOMAGGY_SPECS,
    "legacy_dr8_north": _LEGACY_NANOMAGGY_SPECS,
    "legacy_dr8_south": _LEGACY_NANOMAGGY_SPECS,
}


def magnitude_to_flux_mjy(mag, mag_err, *, zero_jy: float) -> tuple[float, float]:
    """Convert a magnitude and magnitude error to flux density in mJy."""

    mag = _as_float(mag)
    mag_err = _as_float(mag_err)
    if not np.isfinite(mag) or not np.isfinite(mag_err) or mag_err <= 0.0:
        return np.nan, np.nan
    flux_mjy = 1.0e3 * float(zero_jy) * 10.0 ** (-0.4 * mag)
    flux_err_mjy = flux_mjy * np.log(10.0) * 0.4 * mag_err
    return flux_mjy, flux_err_mjy


def nanomaggy_to_flux_mjy(flux_nanomaggy, flux_ivar) -> tuple[float, float]:
    """Convert an AB nanomaggy flux and inverse variance to mJy."""

    flux_nanomaggy = _as_float(flux_nanomaggy)
    flux_ivar = _as_float(flux_ivar)
    if not np.isfinite(flux_nanomaggy) or flux_nanomaggy <= 0.0:
        return np.nan, np.nan
    if not np.isfinite(flux_ivar) or flux_ivar <= 0.0:
        return np.nan, np.nan
    flux_mjy = 0.003631 * flux_nanomaggy
    flux_err_mjy = 0.003631 / np.sqrt(flux_ivar)
    return flux_mjy, flux_err_mjy


def matches_to_photometry(
    matches: Mapping[str, Table],
    *,
    source_id_col: str = "source_id",
    max_mag_err: float | None = None,
    min_quality: int = 2,
) -> Table:
    """Convert raw XMatch outputs into one long-form photometry table.

    Parameters
    ----------
    matches
        Mapping returned by ``bandwagon.xmatch_catalogs``.
    source_id_col
        Identifier column from the uploaded source table.
    max_mag_err
        Optional maximum magnitude uncertainty to keep.
    min_quality
        Minimum numeric flux-quality flag for flux-based catalogs. A value of 2
        keeps confirmed IRAS/AKARI detections while dropping upper limits and
        non-observations.
    """

    tables = [
        _catalog_to_photometry(
            catalog_key,
            table,
            source_id_col=source_id_col,
            max_mag_err=max_mag_err,
            min_quality=min_quality,
        )
        for catalog_key, table in matches.items()
    ]
    tables = [table for table in tables if len(table) > 0]
    if not tables:
        return _empty_photometry_table()
    return vstack(tables, metadata_conflicts="silent")


def _catalog_to_photometry(
    catalog_key: str,
    table: Table,
    *,
    source_id_col: str,
    max_mag_err: float | None,
    min_quality: int,
) -> Table:
    mag_specs = CATALOG_BAND_SPECS.get(catalog_key, ())
    flux_specs = CATALOG_FLUX_SPECS.get(catalog_key, ())
    nanomaggy_specs = CATALOG_NANOMAGGY_SPECS.get(catalog_key, ())
    if not mag_specs and not flux_specs and not nanomaggy_specs:
        return _empty_photometry_table()
    if source_id_col not in table.colnames:
        raise KeyError(f"matched table {catalog_key!r} is missing {source_id_col!r}")

    rows = []
    for row in table:
        source_id = row[source_id_col]
        distance_arcsec = _match_distance_arcsec(row)
        for spec in mag_specs:
            mag = _first_value(row, table.colnames, spec.mag_col)
            mag_err = _first_value(row, table.colnames, spec.err_col)
            if mag is None or mag_err is None:
                continue
            mag = _as_float(mag)
            mag_err = _as_float(mag_err)
            if max_mag_err is not None and np.isfinite(mag_err) and mag_err > max_mag_err:
                continue
            flux_mjy, flux_err_mjy = magnitude_to_flux_mjy(
                mag,
                mag_err,
                zero_jy=spec.zero_jy,
            )
            if not np.isfinite(flux_mjy) or not np.isfinite(flux_err_mjy) or flux_mjy <= 0.0:
                continue
            rows.append(
                {
                    source_id_col: source_id,
                    "catalog": catalog_key,
                    "band": spec.band,
                    "filter_name": spec.filter_name,
                    "speclite_name": SPECLITE_NAMES[spec.filter_name],
                    "mag": mag,
                    "mag_err": mag_err,
                    "mag_system": spec.system,
                    "flux_mjy": flux_mjy,
                    "flux_err_mjy": flux_err_mjy,
                    "psf_fwhm_arcsec": PSF_FWHM_ARCSEC[spec.filter_name],
                    "match_distance_arcsec": distance_arcsec,
                }
            )
        for spec in flux_specs:
            if spec.flux_col not in table.colnames or spec.err_col not in table.colnames:
                continue
            if not _passes_quality(row, spec.quality_col, min_quality=min_quality):
                continue
            flux_jy = _as_float(row[spec.flux_col])
            err_value = _as_float(row[spec.err_col])
            flux_mjy, flux_err_mjy = jy_to_flux_mjy(
                flux_jy,
                err_value,
                err_is_percent=spec.err_is_percent,
            )
            if not np.isfinite(flux_mjy) or not np.isfinite(flux_err_mjy) or flux_mjy <= 0.0:
                continue
            rows.append(
                {
                    source_id_col: source_id,
                    "catalog": catalog_key,
                    "band": spec.band,
                    "filter_name": spec.filter_name,
                    "speclite_name": SPECLITE_NAMES[spec.filter_name],
                    "mag": np.nan,
                    "mag_err": np.nan,
                    "mag_system": "flux",
                    "flux_mjy": flux_mjy,
                    "flux_err_mjy": flux_err_mjy,
                    "psf_fwhm_arcsec": PSF_FWHM_ARCSEC[spec.filter_name],
                    "match_distance_arcsec": distance_arcsec,
                }
            )
        for spec in nanomaggy_specs:
            flux_nanomaggy = _first_value(row, table.colnames, spec.flux_col)
            flux_ivar = _first_value(row, table.colnames, spec.ivar_col)
            if flux_nanomaggy is None or flux_ivar is None:
                continue
            flux_mjy, flux_err_mjy = nanomaggy_to_flux_mjy(flux_nanomaggy, flux_ivar)
            if not np.isfinite(flux_mjy) or not np.isfinite(flux_err_mjy) or flux_mjy <= 0.0:
                continue
            flux_nanomaggy = _as_float(flux_nanomaggy)
            flux_err_nanomaggy = 1.0 / np.sqrt(_as_float(flux_ivar))
            mag = 22.5 - 2.5 * np.log10(flux_nanomaggy)
            mag_err = (2.5 / np.log(10.0)) * (flux_err_nanomaggy / flux_nanomaggy)
            rows.append(
                {
                    source_id_col: source_id,
                    "catalog": catalog_key,
                    "band": spec.band,
                    "filter_name": spec.filter_name,
                    "speclite_name": SPECLITE_NAMES[spec.filter_name],
                    "mag": mag,
                    "mag_err": mag_err,
                    "mag_system": "ab",
                    "flux_mjy": flux_mjy,
                    "flux_err_mjy": flux_err_mjy,
                    "psf_fwhm_arcsec": _row_psf_fwhm(row, table.colnames, spec),
                    "match_distance_arcsec": distance_arcsec,
                }
            )

    if not rows:
        return _empty_photometry_table()
    return _deduplicate_photometry(Table(rows), source_id_col=source_id_col)


def _deduplicate_photometry(table: Table, *, source_id_col: str) -> Table:
    """Keep the best row per source/filter using distance then flux S/N."""

    order = np.lexsort(
        (
            -np.asarray(table["flux_mjy"] / table["flux_err_mjy"], dtype=float),
            np.nan_to_num(np.asarray(table["match_distance_arcsec"], dtype=float), nan=1.0e30),
            np.asarray(table["filter_name"], dtype=str),
            np.asarray(table[source_id_col], dtype=str),
        )
    )
    table = table[order]
    keep = []
    seen = set()
    for i, row in enumerate(table):
        key = (str(row[source_id_col]), str(row["filter_name"]))
        if key in seen:
            continue
        seen.add(key)
        keep.append(i)
    return table[keep]


def _match_distance_arcsec(row) -> float:
    for col in ("angDist", "angDist_arcsec", "_r"):
        if col in row.colnames:
            value = _as_float(row[col])
            return value if np.isfinite(value) else np.nan
    return np.nan


def _first_value(row, colnames, candidates):
    if isinstance(candidates, str):
        candidates = (candidates,)
    for col in candidates:
        if col in colnames:
            return row[col]
    return None


def _row_psf_fwhm(row, colnames, spec: NanomaggyBandSpec) -> float:
    if spec.psf_col is not None:
        value = _first_value(row, colnames, spec.psf_col)
        psf_fwhm = _as_float(value)
        if np.isfinite(psf_fwhm) and psf_fwhm > 0.0:
            return psf_fwhm
    return PSF_FWHM_ARCSEC[spec.filter_name]


def jy_to_flux_mjy(
    flux_jy,
    err_value,
    *,
    err_is_percent: bool = False,
) -> tuple[float, float]:
    """Convert a catalog flux density in Jy to mJy."""

    flux_jy = _as_float(flux_jy)
    err_value = _as_float(err_value)
    if not np.isfinite(flux_jy) or flux_jy <= 0.0:
        return np.nan, np.nan
    if not np.isfinite(err_value) or err_value <= 0.0:
        return np.nan, np.nan
    flux_mjy = 1.0e3 * flux_jy
    if err_is_percent:
        flux_err_mjy = flux_mjy * err_value / 100.0
    else:
        flux_err_mjy = 1.0e3 * err_value
    return flux_mjy, flux_err_mjy


def _passes_quality(row, quality_col: str | None, *, min_quality: int) -> bool:
    if quality_col is None:
        return True
    if quality_col not in row.colnames:
        return True
    quality = _as_float(row[quality_col])
    return np.isfinite(quality) and quality >= min_quality


def _as_float(value) -> float:
    if np.ma.is_masked(value):
        return np.nan
    try:
        return float(value)
    except Exception:
        return np.nan


def _empty_photometry_table() -> Table:
    return Table(
        names=[
            "source_id",
            "catalog",
            "band",
            "filter_name",
            "speclite_name",
            "mag",
            "mag_err",
            "mag_system",
            "flux_mjy",
            "flux_err_mjy",
            "psf_fwhm_arcsec",
            "match_distance_arcsec",
        ],
        dtype=[
            object,
            "U24",
            "U8",
            "U16",
            "U24",
            float,
            float,
            "U8",
            float,
            float,
            float,
            float,
        ],
    )
