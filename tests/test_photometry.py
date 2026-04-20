import numpy as np
from astropy.table import Table

from bandwagon import jy_to_flux_mjy, magnitude_to_flux_mjy, matches_to_photometry


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


def test_matches_to_photometry_converts_default_and_optional_catalogs():
    matches = {
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
        "2mass_xsc": Table(
            {
                "source_id": ["extended"],
                "J.ext": [13.0],
                "e_J.ext": [0.03],
                "H.ext": [12.0],
                "e_H.ext": [0.04],
                "K.ext": [11.0],
                "e_K.ext": [0.05],
                "angDist": [5.0],
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
    }

    phot = matches_to_photometry(matches, min_quality=2)

    assert set(phot["filter_name"]) == {
        "W1",
        "W2",
        "W3",
        "W4",
        "J_2mass",
        "H_2mass",
        "Ks_2mass",
        "S9W_akari",
        "L18W_akari",
        "F12_iras",
    }
    assert np.sum(phot["source_id"] == "extended") == 3
    assert "F25_iras" not in set(phot["filter_name"])
    assert np.all(np.asarray(phot["flux_mjy"], dtype=float) > 0.0)
    assert np.all(np.asarray(phot["flux_err_mjy"], dtype=float) > 0.0)


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
