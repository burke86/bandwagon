import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.providers.sdss import query_sdss_spectra


class FakeSDSS:
    calls = []

    @classmethod
    def query_region(cls, coord, *, spectro, radius, data_release):
        cls.calls.append((coord, spectro, radius, data_release))
        return Table(
            {
                "plate": [1],
                "mjd": [2],
                "fiberID": [3],
                "specObjID": [123],
                "ra": [10.0],
                "dec": [-2.0],
                "z": [0.5],
                "zErr": [0.01],
                "class": ["QSO"],
                "zwarning": [0],
            }
        )


def test_query_sdss_spectra_normalizes_metadata():
    coords = SkyCoord([10.0] * u.deg, [-2.0] * u.deg)
    FakeSDSS.calls = []

    table = query_sdss_spectra(coords, source_id=["src"], client=FakeSDSS, data_release=17)

    assert len(table) == 1
    assert table["source_id"][0] == "src"
    assert table["provider"][0] == "sdss"
    assert table["obs_id"][0] == "1-2-3"
    assert table["redshift"][0] == 0.5
    assert table["spectral_class"][0] == "QSO"
    assert FakeSDSS.calls[0][1] is True
    assert FakeSDSS.calls[0][2] == 2.0 * u.arcsec
    assert FakeSDSS.calls[0][3] == 17
