"""SDSS/BOSS/eBOSS spectra metadata discovery."""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from ..spectra import as_float, as_str, first_col, rows_to_table, source_table_from_coords


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
