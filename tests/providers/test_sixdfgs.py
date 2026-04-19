import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.providers.sixdfgs import query_6dfgs_spectra


def fake_xmatch(coords, catalog, *, radius_arcsec, source_id):
    assert catalog == "VII/259/6dfgs"
    assert radius_arcsec == 5.0
    return Table(
        {
            "source_id": ["src"],
            "Name": ["6dFGS gJ000000-000000"],
            "TargetID": [99],
            "RAJ2000": [10.0],
            "DEJ2000": [-2.0],
            "cz": [30000.0],
            "q_cz": [4],
            "angDist": [1.0],
        }
    )


def test_query_6dfgs_spectra_normalizes_xmatch():
    coords = SkyCoord([10.0] * u.deg, [-2.0] * u.deg)

    table = query_6dfgs_spectra(coords, source_id=["src"], xmatch_func=fake_xmatch)

    assert len(table) == 1
    assert table["provider"][0] == "6dfgs"
    assert table["obs_id"][0] == "99"
    assert table["quality"][0] == "4"
    assert table["redshift"][0] > 0.09
    assert table["sixdfgs_target"][0] == "6dFGS gJ000000-000000"
