"""MAST spectra metadata discovery."""

from __future__ import annotations

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from ..spectra import as_float, as_str, first_col, rows_to_table, source_table_from_coords


def query_mast_spectra(
    coords: SkyCoord,
    *,
    source_id=None,
    radius_arcsec: float = 2.0,
    product_type=("SCIENCE",),
    client=None,
):
    """Query MAST observation products with dataproduct_type='spectrum'."""

    if client is None:
        from astroquery.mast import Observations

        client = Observations

    sources = source_table_from_coords(coords, source_id=source_id)
    rows = []
    for source in sources:
        coord = SkyCoord(source["ra"] * u.deg, source["dec"] * u.deg)
        obs = client.query_region(coord, radius=radius_arcsec * u.arcsec)
        if obs is None or len(obs) == 0:
            continue
        if "dataproduct_type" in obs.colnames:
            obs = obs[np.asarray(obs["dataproduct_type"], dtype=str) == "spectrum"]
        if len(obs) == 0:
            continue
        products = client.get_product_list(obs)
        if products is None or len(products) == 0:
            continue
        if product_type is not None:
            products = client.filter_products(products, productType=list(product_type))
        for row in products:
            rows.append(
                {
                    "source_id": source["source_id"],
                    "provider": "mast",
                    "survey": as_str(first_col(row, ["obs_collection"], "MAST")),
                    "release": "",
                    "instrument": as_str(first_col(row, ["instrument_name"], "")),
                    "obs_id": as_str(first_col(row, ["obs_id", "obsID"], "")),
                    "target_name": as_str(first_col(row, ["target_name"], "")),
                    "ra": as_float(first_col(row, ["s_ra", "ra"], np.nan)),
                    "dec": as_float(first_col(row, ["s_dec", "dec"], np.nan)),
                    "match_distance_arcsec": np.nan,
                    "redshift": np.nan,
                    "redshift_err": np.nan,
                    "spectral_class": "",
                    "quality": as_str(first_col(row, ["calib_level"], "")),
                    "access_url": as_str(first_col(row, ["dataURI", "dataURL"], "")),
                    "access_format": as_str(first_col(row, ["productSubGroupDescription", "extension"], "")),
                    "download_method": "astroquery.mast",
                    "local_path": "",
                    "mast_product_filename": as_str(first_col(row, ["productFilename"], "")),
                }
            )
    return rows_to_table(rows)
