"""LAMOST spectra metadata discovery via VizieR XMatch."""

from __future__ import annotations

import numpy as np
from astropy.coordinates import SkyCoord

from ..spectra import as_float, as_str, first_col, rows_to_table
from ..xmatch import xmatch_catalog


LAMOST_CATALOGS = {
    "lamost_dr5": {
        "catalog": "V/164/dr5",
        "release": "dr5",
    },
}


def query_lamost_spectra(
    coords: SkyCoord,
    *,
    source_id=None,
    radius_arcsec: float = 2.0,
    catalog: str = "lamost_dr5",
    xmatch_func=xmatch_catalog,
):
    """Query LAMOST spectra metadata from the configured VizieR catalog."""

    if catalog not in LAMOST_CATALOGS:
        raise ValueError(f"unknown LAMOST catalog: {catalog!r}")
    config = LAMOST_CATALOGS[catalog]
    matches = xmatch_func(
        coords,
        config["catalog"],
        radius_arcsec=radius_arcsec,
        source_id=source_id,
    )

    rows = []
    for row in matches:
        redshift = as_float(first_col(row, ["z", "Z"], np.nan))
        rows.append(
            {
                "source_id": row["source_id"],
                "provider": "lamost",
                "survey": "LAMOST",
                "release": config["release"],
                "instrument": "LAMOST",
                "obs_id": as_str(first_col(row, ["obsid", "obsID", "ObsID"], "")),
                "target_name": as_str(first_col(row, ["designation", "Designation"], "")),
                "ra": as_float(first_col(row, ["RAJ2000", "ra", "RA"], np.nan)),
                "dec": as_float(first_col(row, ["DEJ2000", "dec", "DEC"], np.nan)),
                "match_distance_arcsec": as_float(first_col(row, ["angDist", "_r"], np.nan)),
                "redshift": redshift if redshift > -1000 else np.nan,
                "redshift_err": as_float(first_col(row, ["z_err", "zerr", "e_z"], np.nan)),
                "spectral_class": as_str(first_col(row, ["class", "Class"], "")),
                "quality": "",
                "access_url": "",
                "access_format": "fits",
                "download_method": "lamost_archive",
                "local_path": "",
                "lamost_obsid": as_str(first_col(row, ["obsid", "obsID", "ObsID"], "")),
            }
        )
    return rows_to_table(rows)
