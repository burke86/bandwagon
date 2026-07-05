"""GAMA spectra metadata discovery via VizieR XMatch."""

from __future__ import annotations

import numpy as np
from astropy.coordinates import SkyCoord

from ..spectra import as_float, as_str, first_col, rows_to_table
from ..xmatch import xmatch_catalog


GAMA_CATALOGS = {
    "gama_dr3": {
        "catalog": "J/MNRAS/474/3875/gamadr3",
        "release": "dr3",
    },
}


def query_gama_spectra(
    coords: SkyCoord,
    *,
    source_id=None,
    radius_arcsec: float = 2.0,
    catalog: str = "gama_dr3",
    xmatch_func=xmatch_catalog,
):
    """Query GAMA spectra metadata from the configured VizieR catalog."""

    if catalog not in GAMA_CATALOGS:
        raise ValueError(f"unknown GAMA catalog: {catalog!r}")
    config = GAMA_CATALOGS[catalog]
    matches = xmatch_func(
        coords,
        config["catalog"],
        radius_arcsec=radius_arcsec,
        source_id=source_id,
    )

    rows = []
    for row in matches:
        redshift = as_float(first_col(row, ["z", "Z"], np.nan))
        quality = as_str(first_col(row, ["q_z", "NQ", "nQ"], ""))
        obs_id = as_str(first_col(row, ["spectID", "SpecId", "SPECID"], ""))
        rows.append(
            {
                "source_id": row["source_id"],
                "provider": "gama",
                "survey": "GAMA",
                "release": config["release"],
                "instrument": "AAOmega",
                "obs_id": obs_id,
                "target_name": as_str(first_col(row, ["CATAID", "Name", "name"], "")),
                "ra": as_float(first_col(row, ["RAJ2000", "ra", "RA"], np.nan)),
                "dec": as_float(first_col(row, ["DEJ2000", "dec", "DEC"], np.nan)),
                "match_distance_arcsec": as_float(first_col(row, ["angDist", "_r"], np.nan)),
                "redshift": redshift if redshift > -1000 else np.nan,
                "redshift_err": as_float(first_col(row, ["e_z", "z_err", "zerr"], np.nan)),
                "spectral_class": as_str(first_col(row, ["class", "Class", "Type"], "")),
                "quality": quality,
                "access_url": "",
                "access_format": "fits",
                "download_method": "gama_archive",
                "local_path": "",
                "gama_spectid": obs_id,
                "gama_survey": as_str(first_col(row, ["Survey", "SURVEY"], "")),
                "gama_survey_code": as_str(first_col(row, ["SCode", "SURVEY_CODE"], "")),
                "gama_wmin_angstrom": as_float(first_col(row, ["Wmin", "WMIN"], np.nan)),
                "gama_wmax_angstrom": as_float(first_col(row, ["Wmax", "WMAX"], np.nan)),
            }
        )
    return rows_to_table(rows)
