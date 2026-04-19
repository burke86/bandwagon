import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon import XMatchError, query_simbad_redshifts, select_best_redshift
from bandwagon.redshift import _simbad_redshift_query


class FakeSimbadClient:
    def __init__(self, tables):
        self.tables = list(tables)
        self.queries = []

    def query_tap(self, query):
        self.queries.append(query)
        return self.tables.pop(0)


def test_simbad_query_contains_redshift_columns_and_cone():
    query = _simbad_redshift_query(
        20.0,
        -10.0,
        radius_arcsec=2.0,
        row_limit=3,
    )

    assert "SELECT TOP 3" in query
    assert "rvz_redshift" in query
    assert "rvz_error" in query
    assert "rvz_type" in query
    assert "rvz_qual" in query
    assert "rvz_bibcode" in query
    assert "CIRCLE('ICRS', 20.000000000000, -10.000000000000" in query


def test_query_simbad_redshifts_normalizes_rows():
    coords = SkyCoord([20.0] * u.deg, [-10.0] * u.deg)
    client = FakeSimbadClient(
        [
            Table(
                {
                    "main_id": ["QSO J0123"],
                    "ra": [20.0001],
                    "dec": [-10.0001],
                    "otype": ["QSO"],
                    "rvz_redshift": [0.047],
                    "rvz_error": [0.001],
                    "rvz_type": ["z"],
                    "rvz_qual": ["A"],
                    "rvz_bibcode": ["2024A&A...1A...1Q"],
                    "match_distance_arcsec": [0.5],
                }
            )
        ]
    )

    redshifts = query_simbad_redshifts(
        coords,
        source_id=["src"],
        radius_arcsec=2.0,
        client=client,
    )

    assert len(redshifts) == 1
    assert redshifts["source_id"][0] == "src"
    assert redshifts["catalog"][0] == "simbad"
    assert redshifts["object_name"][0] == "QSO J0123"
    assert redshifts["redshift"][0] == 0.047
    assert redshifts["redshift_err"][0] == 0.001
    assert redshifts["quality"][0] == "A"
    assert redshifts["match_distance_arcsec"][0] == 0.5
    assert len(client.queries) == 1


def test_query_simbad_redshifts_validates_source_id_length():
    coords = SkyCoord([20.0, 21.0] * u.deg, [-10.0, -11.0] * u.deg)

    with pytest.raises(XMatchError, match="source_id"):
        query_simbad_redshifts(coords, source_id=["src"])


def test_select_best_redshift_prefers_closest_then_quality_then_error():
    redshifts = Table(
        {
            "source_id": ["a", "a", "b", "b"],
            "catalog": ["simbad"] * 4,
            "object_name": ["far_good", "near_ok", "bad_err", "good_err"],
            "redshift": [0.2, 0.21, 0.3, 0.31],
            "redshift_err": [0.001, 0.002, 0.3, 0.01],
            "z_type": ["z"] * 4,
            "quality": ["A", "B", "A", "A"],
            "reference": [""] * 4,
            "object_type": ["QSO"] * 4,
            "match_distance_arcsec": [1.0, 0.2, 0.1, 0.1],
        }
    )

    best = select_best_redshift(redshifts)

    assert best["source_id"].tolist() == ["a", "b"]
    assert best["object_name"].tolist() == ["near_ok", "good_err"]
