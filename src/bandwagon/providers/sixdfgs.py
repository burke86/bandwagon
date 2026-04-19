"""6dFGS spectra metadata discovery via VizieR XMatch."""

from __future__ import annotations

import numpy as np
from astropy.coordinates import SkyCoord

from ..spectra import as_float, as_str, first_col, rows_to_table
from ..xmatch import xmatch_catalog


def query_6dfgs_spectra(
    coords: SkyCoord,
    *,
    source_id=None,
    radius_arcsec: float = 5.0,
    xmatch_func=xmatch_catalog,
):
    """Query 6dFGS final redshift-release metadata from VizieR."""

    matches = xmatch_func(
        coords,
        "VII/259/6dfgs",
        radius_arcsec=radius_arcsec,
        source_id=source_id,
    )

    rows = []
    for row in matches:
        velocity = as_float(first_col(row, ["cz", "HRV", "RadVel", "RV"], np.nan))
        redshift = velocity / 299792.458 if np.isfinite(velocity) else as_float(first_col(row, ["z", "Z"], np.nan))
        rows.append(
            {
                "source_id": row["source_id"],
                "provider": "6dfgs",
                "survey": "6dFGS",
                "release": "final",
                "instrument": "6dF",
                "obs_id": as_str(first_col(row, ["TargetID", "target_id", "ID"], "")),
                "target_name": as_str(first_col(row, ["Name", "name", "6dFGS"], "")),
                "ra": as_float(first_col(row, ["RAJ2000", "ra", "RA"], np.nan)),
                "dec": as_float(first_col(row, ["DEJ2000", "dec", "DEC"], np.nan)),
                "match_distance_arcsec": as_float(first_col(row, ["angDist", "_r"], np.nan)),
                "redshift": redshift,
                "redshift_err": np.nan,
                "spectral_class": as_str(first_col(row, ["SpType", "Type"], "")),
                "quality": as_str(first_col(row, ["q_cz", "Q", "Qual"], "")),
                "access_url": "",
                "access_format": "fits",
                "download_method": "6dfgs_archive",
                "local_path": "",
                "sixdfgs_target": as_str(first_col(row, ["Name", "name", "6dFGS"], "")),
            }
        )
    return rows_to_table(rows)
