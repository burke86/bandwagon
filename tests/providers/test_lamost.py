import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.providers.lamost import query_lamost_spectra


def fake_xmatch(coords, catalog, *, radius_arcsec, source_id):
    assert catalog == "V/164/dr5"
    assert radius_arcsec == 2.0
    return Table(
        {
            "source_id": ["src"],
            "obsid": [42],
            "designation": ["LAMOST J"],
            "ra": [10.0],
            "dec": [-2.0],
            "z": [0.12],
            "z_err": [0.01],
            "class": ["GALAXY"],
            "angDist": [0.4],
        }
    )


def test_query_lamost_spectra_normalizes_xmatch():
    coords = SkyCoord([10.0] * u.deg, [-2.0] * u.deg)

    table = query_lamost_spectra(coords, source_id=["src"], xmatch_func=fake_xmatch)

    assert len(table) == 1
    assert table["provider"][0] == "lamost"
    assert table["release"][0] == "dr5"
    assert table["obs_id"][0] == "42"
    assert table["redshift"][0] == 0.12
    assert table["spectral_class"][0] == "GALAXY"


def test_query_lamost_spectra_rejects_unknown_catalog():
    coords = SkyCoord([10.0] * u.deg, [-2.0] * u.deg)

    with pytest.raises(ValueError, match="unknown LAMOST catalog"):
        query_lamost_spectra(coords, catalog="missing", xmatch_func=fake_xmatch)
