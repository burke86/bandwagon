"""Small wrappers around the CDS XMatch service.

The functions in this module keep the public API narrow: pass an Astropy
``SkyCoord`` array, choose one or more VizieR catalogs, and get matched tables
back. They intentionally do not convert magnitudes to fluxes yet.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astroquery.xmatch import XMatch


SourceInput = SkyCoord | Table


COMMON_CATALOGS: dict[str, str] = {
    "2mass": "II/246/out",
    "2mass_xsc": "VII/233/xsc",
    "akari_fis": "II/298/fis",
    "akari_irc": "II/297/irc",
    "allwise": "II/328/allwise",
    "desi_legacy_dr8_north": "VII/292/north",
    "desi_legacy_dr8_south": "VII/292/south",
    "galex": "II/335/galex_ais",
    "galex_ais": "II/335/galex_ais",
    "iras_psc": "II/125/main",
    "legacy_dr8_north": "VII/292/north",
    "legacy_dr8_south": "VII/292/south",
    "sdss_dr16": "V/154/sdss16",
    "twomass_xsc": "VII/233/xsc",
    "wise": "II/328/allwise",
}

DEFAULT_CATALOGS: dict[str, str] = {
    "galex_ais": "II/335/galex_ais",
    "sdss_dr16": "V/154/sdss16",
    "2mass": "II/246/out",
    "2mass_xsc": "VII/233/xsc",
    "allwise": "II/328/allwise",
    "legacy_dr8_north": "VII/292/north",
    "legacy_dr8_south": "VII/292/south",
}

DEFAULT_RADII_ARCSEC: dict[str, float] = {
    "galex_ais": 3.0,
    "sdss_dr16": 1.0,
    "allwise": 3.0,
    "2mass": 2.0,
    "2mass_xsc": 2.0,
    "legacy_dr8_north": 1.0,
    "legacy_dr8_south": 1.0,
    "akari_irc": 6.0,
    "akari_fis": 20.0,
    "iras_psc": 30.0,
}

CATALOG_BANDS: dict[str, tuple[str, ...]] = {
    "galex_ais": ("FUV", "NUV"),
    "sdss_dr16": ("u", "g", "r", "i", "z"),
    "allwise": ("W1", "W2", "W3", "W4"),
    "2mass": ("J", "H", "Ks"),
    "2mass_xsc": ("J.ext", "H.ext", "K.ext"),
    "akari_irc": ("S9W", "L18W"),
    "akari_fis": ("N60", "WIDE-S", "WIDE-L", "N160"),
    "iras_psc": ("F12", "F25", "F60", "F100"),
    "desi_legacy_dr8_north": ("g", "r", "z", "W1", "W2"),
    "desi_legacy_dr8_south": ("g", "r", "z", "W1", "W2"),
    "legacy_dr8_north": ("photo-z", "value-added"),
    "legacy_dr8_south": ("photo-z", "value-added"),
}


class XMatchError(RuntimeError):
    """Raised when a cross-match request cannot be prepared."""


@dataclass(frozen=True)
class MatchJob:
    """Description of one catalog cross-match."""

    name: str
    catalog: str
    radius_arcsec: float


def make_source_table(
    ra: Iterable[float],
    dec: Iterable[float],
    *,
    source_id: Iterable[object] | None = None,
    id_col: str = "source_id",
    ra_col: str = "ra",
    dec_col: str = "dec",
) -> Table:
    """Build a source table suitable for CDS XMatch.

    Parameters
    ----------
    ra, dec
        ICRS coordinates in degrees.
    source_id
        Optional source identifiers to preserve through the match.
    id_col, ra_col, dec_col
        Output column names.
    """

    ra_values = list(ra)
    dec_values = list(dec)
    if len(ra_values) != len(dec_values):
        raise XMatchError("ra and dec must have the same length.")

    table = Table()
    if source_id is not None:
        source_ids = list(source_id)
        if len(source_ids) != len(ra_values):
            raise XMatchError("source_id must have the same length as coordinates.")
        table[id_col] = source_ids
    table[ra_col] = ra_values
    table[dec_col] = dec_values
    return table


def coords_to_source_table(
    coords: SkyCoord,
    *,
    source_id: Iterable[object] | None = None,
    id_col: str = "source_id",
    ra_col: str = "ra",
    dec_col: str = "dec",
) -> Table:
    """Convert a ``SkyCoord`` scalar or array to a CDS XMatch source table."""

    icrs = coords.icrs
    table = make_source_table(
        np.ravel(np.atleast_1d(icrs.ra.deg)),
        np.ravel(np.atleast_1d(icrs.dec.deg)),
        source_id=source_id,
        id_col=id_col,
        ra_col=ra_col,
        dec_col=dec_col,
    )
    return table


def normalize_vizier_catalog(catalog: str) -> str:
    """Return a CDS XMatch-compatible VizieR catalog identifier."""

    catalog = COMMON_CATALOGS.get(catalog.lower(), catalog)
    if catalog.startswith("vizier:"):
        return catalog
    return f"vizier:{catalog}"


def xmatch_catalog(
    sources: SourceInput,
    catalog: str,
    *,
    radius_arcsec: float,
    source_id: Iterable[object] | None = None,
    id_col: str = "source_id",
    ra_col: str = "ra",
    dec_col: str = "dec",
    cache: bool = True,
) -> Table:
    """Cross-match source coordinates against one VizieR catalog.

    ``sources`` should normally be an Astropy ``SkyCoord`` scalar or array. A
    pre-built ``Table`` with coordinate columns in degrees is also accepted.
    ``catalog`` may be a raw VizieR table id such as ``"II/335/galex_ais"`` or
    a key from ``COMMON_CATALOGS`` such as ``"galex_ais"``.
    """

    source_table = _source_table(
        sources,
        source_id=source_id,
        id_col=id_col,
        ra_col=ra_col,
        dec_col=dec_col,
    )

    if radius_arcsec <= 0:
        raise XMatchError("radius_arcsec must be positive.")

    return XMatch.query(
        cat1=source_table,
        cat2=normalize_vizier_catalog(catalog),
        max_distance=radius_arcsec * u.arcsec,
        colRA1=ra_col,
        colDec1=dec_col,
        cache=cache,
    )


def xmatch_catalogs(
    sources: SourceInput,
    catalogs: Mapping[str, str] | Iterable[str] | None = None,
    *,
    radius_arcsec: float | Mapping[str, float] | None = None,
    source_id: Iterable[object] | None = None,
    id_col: str = "source_id",
    ra_col: str = "ra",
    dec_col: str = "dec",
    cache: bool = True,
) -> dict[str, Table]:
    """Cross-match source coordinates against several VizieR catalogs.

    Parameters
    ----------
    sources
        Astropy ``SkyCoord`` scalar/array, or a pre-built source ``Table``.
    catalogs
        Either a mapping of output name to catalog id, or an iterable of catalog
        ids/common names. Iterable entries are also used as output names. If
        omitted, the default catalogs cover GALEX FUV/NUV, SDSS ugriz, and WISE
        W1/W2/W3/W4.
    radius_arcsec
        One radius for all catalogs, or a mapping of output name to radius. If
        omitted, per-catalog defaults are used.
    """

    if catalogs is None:
        catalogs = DEFAULT_CATALOGS

    if radius_arcsec is None:
        radius_arcsec = DEFAULT_RADII_ARCSEC

    source_table = _source_table(
        sources,
        source_id=source_id,
        id_col=id_col,
        ra_col=ra_col,
        dec_col=dec_col,
    )

    jobs = _match_jobs(catalogs, radius_arcsec=radius_arcsec)
    return {
        job.name: xmatch_catalog(
            source_table,
            job.catalog,
            radius_arcsec=job.radius_arcsec,
            ra_col=ra_col,
            dec_col=dec_col,
            cache=cache,
        )
        for job in jobs
    }


def _source_table(
    sources: SourceInput,
    *,
    source_id: Iterable[object] | None,
    id_col: str,
    ra_col: str,
    dec_col: str,
) -> Table:
    if isinstance(sources, SkyCoord):
        table = coords_to_source_table(
            sources,
            source_id=source_id,
            id_col=id_col,
            ra_col=ra_col,
            dec_col=dec_col,
        )
    elif isinstance(sources, Table):
        if source_id is not None:
            raise XMatchError("source_id can only be passed with SkyCoord input.")
        table = sources
    else:
        raise XMatchError("sources must be an astropy.coordinates.SkyCoord or Table.")

    _validate_sources(table, ra_col=ra_col, dec_col=dec_col)
    return table


def _validate_sources(sources: Table, *, ra_col: str, dec_col: str) -> None:
    missing = [col for col in (ra_col, dec_col) if col not in sources.colnames]
    if missing:
        raise XMatchError(f"source table is missing required columns: {missing}")

    if len(sources) == 0:
        raise XMatchError("source table is empty.")


def _match_jobs(
    catalogs: Mapping[str, str] | Iterable[str],
    *,
    radius_arcsec: float | Mapping[str, float],
) -> list[MatchJob]:
    if isinstance(catalogs, Mapping):
        items = list(catalogs.items())
    else:
        items = [(catalog, catalog) for catalog in catalogs]

    if not items:
        raise XMatchError("at least one catalog is required.")

    jobs: list[MatchJob] = []
    for name, catalog in items:
        radius = _match_radius(name, radius_arcsec)
        if radius <= 0:
            raise XMatchError(f"radius for {name!r} must be positive.")
        jobs.append(MatchJob(name=name, catalog=catalog, radius_arcsec=radius))

    return jobs


def _match_radius(
    name: str,
    radius_arcsec: float | Mapping[str, float],
) -> float:
    if not isinstance(radius_arcsec, Mapping):
        return radius_arcsec

    if name in radius_arcsec:
        return radius_arcsec[name]

    normalized_name = name.lower()
    if normalized_name in radius_arcsec:
        return radius_arcsec[normalized_name]

    raise XMatchError(f"no radius_arcsec value provided for catalog {name!r}.")
