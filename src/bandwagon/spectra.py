"""Archival spectra discovery helpers."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack

from .xmatch import XMatchError


SPECTRA_COLUMNS = [
    "source_id",
    "provider",
    "survey",
    "release",
    "instrument",
    "obs_id",
    "target_name",
    "ra",
    "dec",
    "match_distance_arcsec",
    "redshift",
    "redshift_err",
    "spectral_class",
    "quality",
    "access_url",
    "access_format",
    "download_method",
    "local_path",
]

SPECTRA_DTYPES = [
    object,
    "U24",
    "U32",
    "U24",
    "U32",
    "U64",
    "U96",
    float,
    float,
    float,
    float,
    float,
    "U32",
    "U32",
    "U256",
    "U32",
    "U32",
    "U256",
]


def source_table_from_coords(
    coords: SkyCoord,
    *,
    source_id: Iterable[object] | None = None,
) -> Table:
    """Convert a SkyCoord scalar or array to a source table."""

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


def empty_spectra_table() -> Table:
    """Return an empty normalized spectra-index table."""

    return Table(names=SPECTRA_COLUMNS, dtype=SPECTRA_DTYPES)


def normalize_spectra_tables(tables: Sequence[Table]) -> Table:
    """Stack provider spectra tables using the common schema."""

    tables = [table for table in tables if table is not None and len(table) > 0]
    if not tables:
        return empty_spectra_table()

    normalized = []
    for table in tables:
        data = {}
        n_rows = len(table)
        for col in SPECTRA_COLUMNS:
            if col in table.colnames:
                data[col] = table[col]
            else:
                data[col] = [_default_for_col(col)] * n_rows
        copy = Table(data)
        for col in table.colnames:
            if col not in copy.colnames:
                copy[col] = table[col]
        normalized.append(copy)
    return vstack(normalized, metadata_conflicts="silent")


def query_archival_spectra(
    coords: SkyCoord,
    *,
    source_id: Iterable[object] | None = None,
    providers: Sequence[str] = ("desi", "sdss"),
    radius_arcsec: float = 2.0,
    **provider_kwargs,
) -> Table:
    """Query spectra metadata from one or more providers."""

    if not providers:
        raise ValueError("at least one spectra provider is required.")

    source_table = source_table_from_coords(coords, source_id=source_id)
    source_ids = list(source_table["source_id"])
    results = []
    for provider in providers:
        provider_key = str(provider).lower()
        func = _provider_function(provider_key)
        kwargs = dict(provider_kwargs.get(provider_key, {}))
        results.append(
            func(
                coords,
                source_id=source_ids,
                radius_arcsec=radius_arcsec,
                **kwargs,
            )
        )
    return normalize_spectra_tables(results)


def _provider_function(provider: str):
    if provider == "desi":
        from .providers.desi import query_desi_spectra

        return query_desi_spectra
    if provider == "sdss":
        from .providers.sdss import query_sdss_spectra

        return query_sdss_spectra
    if provider == "lamost":
        from .providers.lamost import query_lamost_spectra

        return query_lamost_spectra
    if provider in {"6dfgs", "sixdfgs"}:
        from .providers.sixdfgs import query_6dfgs_spectra

        return query_6dfgs_spectra
    if provider == "mast":
        from .providers.mast import query_mast_spectra

        return query_mast_spectra
    raise ValueError(f"unknown spectra provider: {provider!r}")


def _default_for_col(col: str):
    if col in {"ra", "dec", "match_distance_arcsec", "redshift", "redshift_err"}:
        return np.nan
    return ""


def as_float(value) -> float:
    """Convert scalar table values to float, returning NaN on failure."""

    if np.ma.is_masked(value):
        return np.nan
    try:
        return float(value)
    except Exception:
        return np.nan


def as_str(value) -> str:
    """Convert scalar table values to stripped strings."""

    if np.ma.is_masked(value):
        return ""
    return str(value).strip()


def first_col(row, names: Sequence[str], default=""):
    """Return the first available row value from a list of possible columns."""

    for name in names:
        if name in row.colnames:
            return row[name]
    return default


def rows_to_table(rows: list[Mapping]) -> Table:
    """Build a table from row dictionaries using the common spectra schema."""

    if not rows:
        return empty_spectra_table()
    return normalize_spectra_tables([Table(rows)])
