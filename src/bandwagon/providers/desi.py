"""DESI spectra metadata discovery via SPARCL."""

from __future__ import annotations

import numpy as np
from astropy.coordinates import SkyCoord

from ..spectra import as_float, as_str, rows_to_table, source_table_from_coords


DESI_OUTFIELDS = [
    "sparcl_id",
    "targetid",
    "ra",
    "dec",
    "redshift",
    "redshift_err",
    "redshift_warning",
    "spectype",
    "survey",
    "instrument",
    "data_release",
    "wavemin",
    "wavemax",
    "exptime",
    "specprimary",
]


def query_desi_spectra(
    coords: SkyCoord,
    *,
    source_id=None,
    radius_arcsec: float = 2.0,
    release: str = "dr1",
    row_limit: int = 5,
    client=None,
):
    """Query DESI spectra metadata from SPARCL.

    This provider uses a small RA/Dec bounding box query and filters matches by
    angular separation locally.
    """

    if client is None:
        try:
            from sparcl.client import SparclClient
        except ImportError as exc:
            raise ImportError("Install bandwagon[spectra] to query DESI/SPARCL spectra.") from exc

        client = SparclClient()

    sources = source_table_from_coords(coords, source_id=source_id)
    rows = []
    for source in sources:
        result = _find_desi_records(
            client,
            float(source["ra"]),
            float(source["dec"]),
            radius_arcsec=radius_arcsec,
            release=release,
            row_limit=row_limit,
        )
        records = getattr(result, "records", result)
        if records is None:
            continue
        for record in records:
            distance = _distance_arcsec(source["ra"], source["dec"], record)
            if not np.isfinite(distance) or distance > radius_arcsec:
                continue
            rows.append(
                {
                    "source_id": source["source_id"],
                    "provider": "desi",
                    "survey": as_str(record.get("survey", "DESI")),
                    "release": as_str(record.get("data_release", release)),
                    "instrument": as_str(record.get("instrument", "DESI")),
                    "obs_id": as_str(record.get("sparcl_id", "")),
                    "target_name": as_str(record.get("targetid", "")),
                    "ra": as_float(record.get("ra", np.nan)),
                    "dec": as_float(record.get("dec", np.nan)),
                    "match_distance_arcsec": distance,
                    "redshift": as_float(record.get("redshift", np.nan)),
                    "redshift_err": as_float(record.get("redshift_err", np.nan)),
                    "spectral_class": as_str(record.get("spectype", "")),
                    "quality": as_str(record.get("redshift_warning", "")),
                    "access_url": "",
                    "access_format": "sparcl",
                    "download_method": "sparcl",
                    "local_path": "",
                    "desi_targetid": as_str(record.get("targetid", "")),
                    "desi_sparcl_id": as_str(record.get("sparcl_id", "")),
                    "desi_redshift_warning": as_str(record.get("redshift_warning", "")),
                }
            )
    return rows_to_table(rows)


def _find_desi_records(
    client,
    ra_deg: float,
    dec_deg: float,
    *,
    radius_arcsec: float,
    release: str,
    row_limit: int,
):
    radius_deg = radius_arcsec / 3600.0
    constraints = {
        "ra": [ra_deg - radius_deg, ra_deg + radius_deg],
        "dec": [dec_deg - radius_deg, dec_deg + radius_deg],
    }
    if release:
        constraints["data_release"] = [release]
    return client.find(
        outfields=DESI_OUTFIELDS,
        constraints=constraints,
        limit=row_limit,
    )


def _distance_arcsec(ra_deg, dec_deg, record) -> float:
    ra = as_float(record.get("ra", np.nan))
    dec = as_float(record.get("dec", np.nan))
    if not np.isfinite(ra) or not np.isfinite(dec):
        return np.nan
    c0 = SkyCoord(ra_deg, dec_deg, unit="deg")
    c1 = SkyCoord(ra, dec, unit="deg")
    return float(c0.separation(c1).arcsec)
