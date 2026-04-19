import astropy.units as u
from astropy.coordinates import SkyCoord

from bandwagon.providers.desi import query_desi_spectra


class Result:
    def __init__(self, records):
        self.records = records


class FakeSparclClient:
    def __init__(self):
        self.calls = []

    def find(self, *, outfields, constraints, limit):
        self.calls.append((outfields, constraints, limit))
        return Result(
            [
                {
                    "sparcl_id": "abc",
                    "targetid": 123,
                    "ra": 10.0001,
                    "dec": -2.0,
                    "redshift": 1.2,
                    "redshift_err": 0.01,
                    "redshift_warning": 0,
                    "spectype": "QSO",
                    "survey": "main",
                    "instrument": "DESI",
                    "data_release": "dr1",
                }
            ]
        )


def test_query_desi_spectra_normalizes_sparcl_records():
    coords = SkyCoord([10.0] * u.deg, [-2.0] * u.deg)
    client = FakeSparclClient()

    table = query_desi_spectra(coords, source_id=["src"], client=client, radius_arcsec=2.0)

    assert len(table) == 1
    assert table["provider"][0] == "desi"
    assert table["obs_id"][0] == "abc"
    assert table["target_name"][0] == "123"
    assert table["redshift"][0] == 1.2
    assert table["spectral_class"][0] == "QSO"
    assert table["match_distance_arcsec"][0] < 1.0
    assert client.calls[0][2] == 5
