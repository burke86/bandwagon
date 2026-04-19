"""SIMBAD redshift helpers."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astroquery.simbad import Simbad

from .xmatch import XMatchError


def query_simbad_redshifts(
    coords: SkyCoord,
    *,
    source_id: Iterable[object] | None = None,
    radius_arcsec: float = 2.0,
    row_limit: int = 5,
    client=None,
) -> Table:
    """Query SIMBAD TAP for redshifts near each coordinate.

    The result is a long-form table with one or more redshift candidates per
    input coordinate. Use ``select_best_redshift`` to keep one row per source.
    """

    if radius_arcsec <= 0:
        raise XMatchError("radius_arcsec must be positive.")
    if row_limit <= 0:
        raise XMatchError("row_limit must be positive.")

    source_table = _source_table_from_coords(coords, source_id=source_id)
    simbad = Simbad if client is None else client
    tables = []
    for row in source_table:
        query = _simbad_redshift_query(
            float(row["ra"]),
            float(row["dec"]),
            radius_arcsec=radius_arcsec,
            row_limit=row_limit,
        )
        result = simbad.query_tap(query)
        if result is None or len(result) == 0:
            continue
        tables.append(
            _normalize_simbad_rows(
                result,
                source_id=row["source_id"],
            )
        )

    if not tables:
        return _empty_redshift_table()
    return vstack(tables, metadata_conflicts="silent")


def select_best_redshift(redshifts: Table) -> Table:
    """Select the preferred SIMBAD redshift candidate for each source."""

    if len(redshifts) == 0:
        return _empty_redshift_table()

    order = np.lexsort(
        (
            np.nan_to_num(np.asarray(redshifts["redshift_err"], dtype=float), nan=1.0e30),
            np.asarray([_quality_rank(q) for q in redshifts["quality"]], dtype=int),
            np.nan_to_num(np.asarray(redshifts["match_distance_arcsec"], dtype=float), nan=1.0e30),
            np.asarray(redshifts["source_id"], dtype=str),
        )
    )
    sorted_redshifts = redshifts[order]

    keep = []
    seen = set()
    for i, row in enumerate(sorted_redshifts):
        source = str(row["source_id"])
        if source in seen:
            continue
        seen.add(source)
        keep.append(i)
    return sorted_redshifts[keep]


def _source_table_from_coords(
    coords: SkyCoord,
    *,
    source_id: Iterable[object] | None,
) -> Table:
    if not isinstance(coords, SkyCoord):
        raise XMatchError("coords must be an astropy.coordinates.SkyCoord.")

    icrs = coords.icrs
    ra = np.ravel(np.atleast_1d(icrs.ra.deg))
    dec = np.ravel(np.atleast_1d(icrs.dec.deg))
    if source_id is None:
        source_ids = [f"source_{i}" for i in range(len(ra))]
    else:
        source_ids = list(source_id)
        if len(source_ids) != len(ra):
            raise XMatchError("source_id must have the same length as coords.")

    return Table({"source_id": source_ids, "ra": ra, "dec": dec})


def _simbad_redshift_query(
    ra_deg: float,
    dec_deg: float,
    *,
    radius_arcsec: float,
    row_limit: int,
) -> str:
    radius_deg = radius_arcsec / 3600.0
    return f"""
SELECT TOP {int(row_limit)}
    main_id,
    ra,
    dec,
    otype,
    rvz_redshift,
    rvz_error,
    rvz_type,
    rvz_qual,
    rvz_bibcode,
    DISTANCE(
        POINT('ICRS', basic.ra, basic.dec),
        POINT('ICRS', {ra_deg:.12f}, {dec_deg:.12f})
    ) * 3600.0 AS match_distance_arcsec
FROM basic
WHERE rvz_redshift IS NOT NULL
AND CONTAINS(
    POINT('ICRS', basic.ra, basic.dec),
    CIRCLE('ICRS', {ra_deg:.12f}, {dec_deg:.12f}, {radius_deg:.12f})
) = 1
ORDER BY match_distance_arcsec ASC
"""


def _normalize_simbad_rows(rows: Table, *, source_id) -> Table:
    out = []
    for row in rows:
        redshift = _as_float(row["rvz_redshift"])
        if not np.isfinite(redshift):
            continue
        out.append(
            {
                "source_id": source_id,
                "catalog": "simbad",
                "object_name": _as_str(row, "main_id"),
                "redshift": redshift,
                "redshift_err": _as_float(row["rvz_error"]) if "rvz_error" in row.colnames else np.nan,
                "z_type": _as_str(row, "rvz_type"),
                "quality": _as_str(row, "rvz_qual"),
                "reference": _as_str(row, "rvz_bibcode"),
                "object_type": _as_str(row, "otype"),
                "match_distance_arcsec": _as_float(row["match_distance_arcsec"])
                if "match_distance_arcsec" in row.colnames
                else np.nan,
            }
        )
    if not out:
        return _empty_redshift_table()
    return Table(out)


def _quality_rank(value) -> int:
    text = str(value).strip().upper()
    if len(text) == 1 and "A" <= text <= "E":
        return ord(text) - ord("A")
    return 99


def _as_float(value) -> float:
    if np.ma.is_masked(value):
        return np.nan
    try:
        return float(value)
    except Exception:
        return np.nan


def _as_str(row, col: str) -> str:
    if col not in row.colnames:
        return ""
    value = row[col]
    if np.ma.is_masked(value):
        return ""
    return str(value).strip()


def _empty_redshift_table() -> Table:
    return Table(
        names=[
            "source_id",
            "catalog",
            "object_name",
            "redshift",
            "redshift_err",
            "z_type",
            "quality",
            "reference",
            "object_type",
            "match_distance_arcsec",
        ],
        dtype=[
            object,
            "U16",
            "U64",
            float,
            float,
            "U8",
            "U8",
            "U32",
            "U32",
            float,
        ],
    )
