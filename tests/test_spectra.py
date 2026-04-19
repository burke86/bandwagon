import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.spectra import normalize_spectra_tables, query_archival_spectra, source_table_from_coords


def test_source_table_from_coords_validates_source_ids():
    coords = SkyCoord([1, 2] * u.deg, [3, 4] * u.deg)

    with pytest.raises(Exception, match="source_id"):
        source_table_from_coords(coords, source_id=["one"])


def test_normalize_spectra_tables_preserves_common_and_extra_columns():
    table = Table(
        {
            "source_id": ["a"],
            "provider": ["test"],
            "redshift": [0.1],
            "extra_col": ["extra"],
        }
    )

    out = normalize_spectra_tables([table])

    assert out["source_id"][0] == "a"
    assert out["provider"][0] == "test"
    assert out["redshift"][0] == 0.1
    assert out["extra_col"][0] == "extra"
    assert "access_url" in out.colnames


def test_query_archival_spectra_rejects_unknown_provider():
    coords = SkyCoord([1] * u.deg, [3] * u.deg)

    with pytest.raises(ValueError, match="unknown spectra provider"):
        query_archival_spectra(coords, providers=("missing",))
