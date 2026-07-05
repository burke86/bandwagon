import numpy as np
from astropy.table import Table

from bandwagon import (
    jy_to_flux_mjy,
    magnitude_to_flux_mjy,
    matches_to_photometry,
    nanomaggy_to_flux_mjy,
)


def test_magnitude_to_flux_mjy_uses_supplied_zero_point():
    flux, err = magnitude_to_flux_mjy(0.0, 0.1, zero_jy=3631.0)

    assert flux == 3_631_000.0
    assert np.isclose(err, flux * np.log(10.0) * 0.4 * 0.1)


def test_jy_to_flux_mjy_handles_absolute_and_percent_errors():
    flux, err = jy_to_flux_mjy(0.5, 0.05)
    assert flux == 500.0
    assert err == 50.0

    flux, err = jy_to_flux_mjy(0.5, 10.0, err_is_percent=True)
    assert flux == 500.0
    assert err == 50.0


def test_nanomaggy_to_flux_mjy_uses_inverse_variance():
    flux, err = nanomaggy_to_flux_mjy(1000.0, 25.0)

    assert np.isclose(flux, 3.631)
    assert np.isclose(err, 0.003631 / 5.0)

    flux, err = nanomaggy_to_flux_mjy(-1.0, 25.0)
    assert np.isnan(flux)
    assert np.isnan(err)


def test_matches_to_photometry_converts_default_and_optional_catalogs():
    matches = {
        "sdss_dr16": Table(
            {
                "source_id": ["src"],
                "umag": [10.0],
                "e_umag": [0.01],
                "gmag": [10.0],
                "e_gmag": [0.01],
                "rmag": [10.0],
                "e_rmag": [0.01],
                "imag": [10.0],
                "e_imag": [0.01],
                "zmag": [10.0],
                "e_zmag": [0.01],
                "upmag": [16.0],
                "e_upmag": [0.05],
                "gpmag": [15.5],
                "e_gpmag": [0.05],
                "rpmag": [15.0],
                "e_rpmag": [0.05],
                "ipmag": [14.8],
                "e_ipmag": [0.05],
                "zpmag": [14.7],
                "e_zpmag": [0.05],
                "angDist": [0.1],
            }
        ),
        "allwise": Table(
            {
                "source_id": ["src"],
                "W1mag": [15.0],
                "e_W1mag": [0.04],
                "W2mag": [14.8],
                "e_W2mag": [0.05],
                "W3mag": [12.0],
                "e_W3mag": [0.2],
                "W4mag": [9.0],
                "e_W4mag": [0.3],
                "angDist": [0.3],
            }
        ),
        "2mass": Table(
            {
                "source_id": ["src"],
                "Jmag": [16.0],
                "e_Jmag": [0.05],
                "Hmag": [15.5],
                "e_Hmag": [0.06],
                "Kmag": [15.0],
                "e_Kmag": [0.07],
                "angDist": [0.2],
            }
        ),
        "akari_irc": Table(
            {
                "source_id": ["src"],
                "S09": [0.1],
                "e_S09": [0.01],
                "q_S09": [3],
                "S18": [0.2],
                "e_S18": [0.02],
                "q_S18": [3],
            }
        ),
        "iras_psc": Table(
            {
                "source_id": ["src"],
                "Fnu_12": [0.5],
                "e_Fnu_12": [10.0],
                "q_Fnu_12": [3],
                "Fnu_25": [0.6],
                "e_Fnu_25": [20.0],
                "q_Fnu_25": [1],
            }
        ),
        "legacy_dr8_north": Table(
            {
                "source_id": ["src"],
                "FLUX_G": [1000.0],
                "FLUX_IVAR_G": [25.0],
                "FLUX_R": [800.0],
                "FLUX_IVAR_R": [16.0],
                "FLUX_Z": [-1.0],
                "FLUX_IVAR_Z": [9.0],
                "FLUX_W1": [500.0],
                "FLUX_IVAR_W1": [4.0],
                "FLUX_W2": [400.0],
                "FLUX_IVAR_W2": [1.0],
                "PSFSIZE_G": [1.4],
                "PSFSIZE_R": [1.3],
                "angDist": [0.15],
            }
        ),
    }

    phot = matches_to_photometry(matches, min_quality=2)

    assert set(phot["filter_name"]) == {
        "W1",
        "W2",
        "W3",
        "W4",
        "u_sdss",
        "g_sdss",
        "r_sdss",
        "i_sdss",
        "z_sdss",
        "J_2mass",
        "H_2mass",
        "Ks_2mass",
        "S9W_akari",
        "L18W_akari",
        "F12_iras",
        "g_legacy",
        "r_legacy",
        "W1_legacy",
        "W2_legacy",
    }
    assert "z_legacy" not in set(phot["filter_name"])
    assert "F25_iras" not in set(phot["filter_name"])
    assert np.all(np.asarray(phot["flux_mjy"], dtype=float) > 0.0)
    assert np.all(np.asarray(phot["flux_err_mjy"], dtype=float) > 0.0)

    sdss_u = phot[phot["filter_name"] == "u_sdss"][0]
    assert np.isclose(sdss_u["mag"], 16.0)
    assert np.isclose(sdss_u["flux_mjy"], 1.0e3 * 3631.0 * 10.0 ** (-0.4 * 16.0))
    assert sdss_u["photometry_method"] == "psf"

    legacy_g = phot[phot["filter_name"] == "g_legacy"][0]
    assert np.isclose(legacy_g["flux_mjy"], 3.631)
    assert np.isclose(legacy_g["mag"], 15.0)
    assert np.isclose(legacy_g["psf_fwhm_arcsec"], 1.4)
    assert legacy_g["photometry_method"] == "catalog"


def test_matches_to_photometry_deduplicates_by_distance_then_snr():
    matches = {
        "akari_irc": Table(
            {
                "source_id": ["src", "src"],
                "S09": [0.1, 0.2],
                "e_S09": [0.01, 0.01],
                "q_S09": [3, 3],
                "angDist": [2.0, 1.0],
            }
        )
    }

    phot = matches_to_photometry(matches)

    assert len(phot) == 1
    assert phot["filter_name"][0] == "S9W_akari"
    assert phot["flux_mjy"][0] == 200.0
