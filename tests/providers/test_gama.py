import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.providers.gama import query_gama_spectra


def fake_xmatch(coords, catalog, *, radius_arcsec, source_id):
    assert catalog == "J/MNRAS/474/3875/gamadr3"
    assert radius_arcsec == 2.0
    return Table(
        {
            "source_id": ["src"],
            "spectID": ["G12_Y1_001_001"],
            "Survey": ["GAMA"],
            "SCode": [1],
            "RAJ2000": [180.0],
            "DEJ2000": [0.5],
            "Wmin": [3720.0],
            "Wmax": [8850.0],
            "z": [0.22],
            "q_z": [4],
            "angDist": [0.6],
        }
    )


def test_query_gama_spectra_normalizes_xmatch():
    coords = SkyCoord([180.0] * u.deg, [0.5] * u.deg)

    table = query_gama_spectra(coords, source_id=["src"], xmatch_func=fake_xmatch)

    assert len(table) == 1
    assert table["provider"][0] == "gama"
    assert table["release"][0] == "dr3"
    assert table["instrument"][0] == "AAOmega"
    assert table["obs_id"][0] == "G12_Y1_001_001"
    assert table["redshift"][0] == 0.22
    assert table["quality"][0] == "4"
    assert table["gama_wmin_angstrom"][0] == 3720.0
    assert table["gama_wmax_angstrom"][0] == 8850.0


def test_query_gama_spectra_rejects_unknown_catalog():
    coords = SkyCoord([180.0] * u.deg, [0.5] * u.deg)

    with pytest.raises(ValueError, match="unknown GAMA catalog"):
        query_gama_spectra(coords, catalog="missing", xmatch_func=fake_xmatch)
