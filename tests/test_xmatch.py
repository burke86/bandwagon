import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon import (
    DEFAULT_CATALOGS,
    XMatchError,
    coords_to_source_table,
    normalize_vizier_catalog,
    xmatch_catalogs,
)
import bandwagon.xmatch as xmatch_module
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


def test_default_catalogs_include_2mass():
    assert DEFAULT_CATALOGS["2mass"] == "II/246/out"
    assert DEFAULT_CATALOGS["legacy_dr8_north"] == "VII/292/north"
    assert DEFAULT_CATALOGS["legacy_dr8_south"] == "VII/292/south"


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


def test_xmatch_catalogs_enriches_sdss_psf_by_default(monkeypatch):
    calls = []

    def fake_xmatch_catalog(sources, catalog, *, radius_arcsec, ra_col, dec_col, cache):
        calls.append((catalog, radius_arcsec))
        return Table({"source_id": ["src"], "objID": [123], "ra": [10.0], "dec": [-2.0]})

    def fake_enrich(table, *, chunk_size, cache):
        out = table.copy()
        out["upmag"] = [18.0]
        out["e_upmag"] = [0.01]
        out["chunk_size"] = [chunk_size]
        return out

    monkeypatch.setattr(xmatch_module, "xmatch_catalog", fake_xmatch_catalog)
    import bandwagon.providers.sdss as sdss_provider

    monkeypatch.setattr(sdss_provider, "enrich_sdss_psf_photometry", fake_enrich)

    result = xmatch_catalogs(
        SkyCoord([10.0] * u.deg, [-2.0] * u.deg),
        catalogs={"sdss_dr16": "V/154/sdss16", "2mass": "II/246/out"},
        radius_arcsec={"sdss_dr16": 1.0, "2mass": 2.0},
        source_id=["src"],
        sdss_psf_chunk_size=25,
    )

    assert calls == [("V/154/sdss16", 1.0), ("II/246/out", 2.0)]
    assert result["sdss_dr16"]["upmag"][0] == 18.0
    assert result["sdss_dr16"]["chunk_size"][0] == 25
    assert "upmag" not in result["2mass"].colnames
