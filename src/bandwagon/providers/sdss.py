"""SDSS/BOSS/eBOSS metadata and PSF photometry discovery."""

from __future__ import annotations

from collections.abc import Sequence

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack

from ..spectra import as_float, as_str, first_col, rows_to_table, source_table_from_coords

SDSS_PSF_COLUMNS = (
    "upmag",
    "gpmag",
    "rpmag",
    "ipmag",
    "zpmag",
    "e_upmag",
    "e_gpmag",
    "e_rpmag",
    "e_ipmag",
    "e_zpmag",
)


def enrich_sdss_psf_photometry(
    table: Table,
    *,
    catalog: str = "V/154/sdss16",
    chunk_size: int = 500,
    cache: bool = True,
    vizier_cls=None,
) -> Table:
    """Attach SDSS PSF magnitude columns to an XMatch table using ``objID``.

    CDS XMatch exposes only the compact SDSS model-magnitude view for
    ``V/154/sdss16``. The full VizieR table has PSF columns, so we use the
    XMatch ``objID`` results as a fast spatial pre-match and fetch PSF columns
    by identifier in batches.
    """

    out = table.copy(copy_data=True)
    if len(out) == 0 or "objID" not in out.colnames:
        return out
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive.")

    objids = _unique_objids(out["objID"])
    if not objids:
        return out

    psf_rows = _query_sdss_psf_by_objid(
        objids,
        catalog=catalog,
        chunk_size=chunk_size,
        cache=cache,
        vizier_cls=vizier_cls,
    )
    if len(psf_rows) == 0:
        return out

    psf_by_objid = {_objid_key(row["objID"]): row for row in psf_rows if _objid_key(row["objID"])}
    for col in SDSS_PSF_COLUMNS:
        values = np.full(len(out), np.nan, dtype=float)
        for i, row in enumerate(out):
            psf_row = psf_by_objid.get(_objid_key(row["objID"]))
            if psf_row is not None and col in psf_row.colnames:
                values[i] = as_float(psf_row[col])
        if col in out.colnames:
            out[col] = values
        else:
            out.add_column(values, name=col)
    return out


def _query_sdss_psf_by_objid(
    objids: Sequence[str],
    *,
    catalog: str,
    chunk_size: int,
    cache: bool,
    vizier_cls=None,
) -> Table:
    if vizier_cls is None:
        from astroquery.vizier import Vizier

        vizier_cls = Vizier

    tables = []
    columns = ["objID", *SDSS_PSF_COLUMNS]
    for chunk in _chunks(list(objids), chunk_size):
        vizier = vizier_cls(columns=columns, row_limit=-1)
        result = vizier.query_constraints(catalog=catalog, objID=",".join(chunk), cache=cache)
        if result:
            tables.append(result[0])
    if not tables:
        return Table(names=columns, dtype=[object, *([float] * len(SDSS_PSF_COLUMNS))])
    return vstack(tables, metadata_conflicts="silent")


def _unique_objids(values) -> list[str]:
    seen = set()
    objids = []
    for value in values:
        key = _objid_key(value)
        if key and key not in seen:
            seen.add(key)
            objids.append(key)
    return objids


def _objid_key(value) -> str:
    if np.ma.is_masked(value):
        return ""
    text = str(value).strip()
    if not text or text == "--":
        return ""
    try:
        return str(int(text))
    except ValueError:
        return text


def _chunks(values: list[str], chunk_size: int):
    for start in range(0, len(values), chunk_size):
        yield values[start : start + chunk_size]


def query_sdss_spectra(
    coords: SkyCoord,
    *,
    source_id=None,
    radius_arcsec: float = 2.0,
    data_release: int = 17,
    client=None,
):
    """Query SDSS spectra metadata near each coordinate."""

    if client is None:
        from astroquery.sdss import SDSS

        client = SDSS

    sources = source_table_from_coords(coords, source_id=source_id)
    rows = []
    for source in sources:
        coord = SkyCoord(source["ra"] * u.deg, source["dec"] * u.deg)
        result = client.query_region(
            coord,
            spectro=True,
            radius=radius_arcsec * u.arcsec,
            data_release=data_release,
        )
        if result is None:
            continue
        for row in result:
            redshift = as_float(first_col(row, ["z", "Z"], np.nan))
            rows.append(
                {
                    "source_id": source["source_id"],
                    "provider": "sdss",
                    "survey": "SDSS",
                    "release": f"dr{data_release}",
                    "instrument": "SDSS",
                    "obs_id": _sdss_obs_id(row),
                    "target_name": as_str(first_col(row, ["specObjID", "specobjid"], "")),
                    "ra": as_float(first_col(row, ["ra", "RA"], np.nan)),
                    "dec": as_float(first_col(row, ["dec", "DEC"], np.nan)),
                    "match_distance_arcsec": np.nan,
                    "redshift": redshift if redshift > -1000 else np.nan,
                    "redshift_err": as_float(first_col(row, ["zErr", "zerr", "Z_ERR"], np.nan)),
                    "spectral_class": as_str(first_col(row, ["class", "CLASS"], "")),
                    "quality": as_str(first_col(row, ["zwarning", "ZWARNING"], "")),
                    "access_url": "",
                    "access_format": "fits",
                    "download_method": "astroquery.sdss",
                    "local_path": "",
                    "sdss_plate": as_str(first_col(row, ["plate", "PLATE"], "")),
                    "sdss_mjd": as_str(first_col(row, ["mjd", "MJD"], "")),
                    "sdss_fiberid": as_str(first_col(row, ["fiberID", "fiberid", "FIBERID"], "")),
                    "sdss_specobjid": as_str(first_col(row, ["specObjID", "specobjid"], "")),
                }
            )
    return rows_to_table(rows)


def _sdss_obs_id(row) -> str:
    plate = as_str(first_col(row, ["plate", "PLATE"], ""))
    mjd = as_str(first_col(row, ["mjd", "MJD"], ""))
    fiber = as_str(first_col(row, ["fiberID", "fiberid", "FIBERID"], ""))
    if plate and mjd and fiber:
        return f"{plate}-{mjd}-{fiber}"
    return as_str(first_col(row, ["specObjID", "specobjid"], ""))
