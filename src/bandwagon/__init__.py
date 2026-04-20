"""Bulk photometric catalog cross-matching helpers."""

from .photometry import (
    CATALOG_BAND_SPECS,
    CATALOG_FLUX_SPECS,
    CATALOG_NANOMAGGY_SPECS,
    PSF_FWHM_ARCSEC,
    SPECLITE_NAMES,
    TWOMASS_VEGA_ZEROPOINT_JY,
    WISE_VEGA_ZEROPOINT_JY,
    jy_to_flux_mjy,
    magnitude_to_flux_mjy,
    matches_to_photometry,
    nanomaggy_to_flux_mjy,
)
from .redshift import query_simbad_redshifts, select_best_redshift
from .spectra import empty_spectra_table, query_archival_spectra
from .xmatch import (
    CATALOG_BANDS,
    COMMON_CATALOGS,
    DEFAULT_CATALOGS,
    DEFAULT_RADII_ARCSEC,
    XMatchError,
    coords_to_source_table,
    make_source_table,
    normalize_vizier_catalog,
    xmatch_catalog,
    xmatch_catalogs,
)

__all__ = [
    "CATALOG_BAND_SPECS",
    "CATALOG_BANDS",
    "CATALOG_FLUX_SPECS",
    "CATALOG_NANOMAGGY_SPECS",
    "COMMON_CATALOGS",
    "DEFAULT_CATALOGS",
    "DEFAULT_RADII_ARCSEC",
    "PSF_FWHM_ARCSEC",
    "SPECLITE_NAMES",
    "TWOMASS_VEGA_ZEROPOINT_JY",
    "WISE_VEGA_ZEROPOINT_JY",
    "XMatchError",
    "coords_to_source_table",
    "empty_spectra_table",
    "jy_to_flux_mjy",
    "magnitude_to_flux_mjy",
    "make_source_table",
    "matches_to_photometry",
    "nanomaggy_to_flux_mjy",
    "normalize_vizier_catalog",
    "query_simbad_redshifts",
    "query_archival_spectra",
    "select_best_redshift",
    "xmatch_catalog",
    "xmatch_catalogs",
]
