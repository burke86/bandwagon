import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon import (
    XMatchError,
    coords_to_source_table,
    normalize_vizier_catalog,
    xmatch_catalogs,
)
from bandwagon.xmatch import _match_jobs


def test_coords_to_source_table_uses_skycoord_icrs_degrees():
    coords = SkyCoord([10.0, 11.0] * u.deg, [-2.0, -3.0] * u.deg, frame="icrs")

    table = coords_to_source_table(coords, source_id=["a", "b"])

    assert table.colnames == ["source_id", "ra", "dec"]
    assert table["source_id"].tolist() == ["a", "b"]
    assert table["ra"].tolist() == [10.0, 11.0]
    assert table["dec"].tolist() == [-2.0, -3.0]


def test_coords_to_source_table_validates_source_id_length():
    coords = SkyCoord([10.0, 11.0] * u.deg, [-2.0, -3.0] * u.deg, frame="icrs")

    with pytest.raises(XMatchError, match="source_id"):
        coords_to_source_table(coords, source_id=["a"])


def test_normalize_vizier_catalog_resolves_aliases():
    assert normalize_vizier_catalog("2mass") == "vizier:II/246/out"
    assert normalize_vizier_catalog("akari_irc") == "vizier:II/297/irc"
    assert normalize_vizier_catalog("akari_fis") == "vizier:II/298/fis"
    assert normalize_vizier_catalog("iras_psc") == "vizier:II/125/main"
    assert normalize_vizier_catalog("vizier:II/246/out") == "vizier:II/246/out"


def test_match_jobs_use_default_optional_radii():
    jobs = _match_jobs(["2mass", "akari_irc", "akari_fis", "iras_psc"], radius_arcsec={
        "2mass": 2.0,
        "akari_irc": 6.0,
        "akari_fis": 20.0,
        "iras_psc": 30.0,
    })

    assert [(job.name, job.radius_arcsec) for job in jobs] == [
        ("2mass", 2.0),
        ("akari_irc", 6.0),
        ("akari_fis", 20.0),
        ("iras_psc", 30.0),
    ]


def test_xmatch_catalogs_rejects_source_id_with_table_input():
    table = Table({"ra": [10.0], "dec": [-2.0]})

    with pytest.raises(XMatchError, match="source_id"):
        xmatch_catalogs(table, source_id=["a"])
