import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.providers.sdss import enrich_sdss_psf_photometry, query_sdss_spectra


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


class FakeVizier:
    calls = []

    def __init__(self, *, columns, row_limit):
        self.columns = columns
        self.row_limit = row_limit

    def query_constraints(self, *, catalog, objID, cache):
        self.__class__.calls.append((catalog, objID, cache, tuple(self.columns), self.row_limit))
        rows = {
            "123": (18.0, 17.0, 16.0, 15.0, 14.0),
            "456": (19.0, 18.0, 17.0, 16.0, 15.0),
        }
        selected = [key for key in objID.split(",") if key in rows]
        if not selected:
            return []
        return [
            Table(
                {
                    "objID": selected,
                    "upmag": [rows[key][0] for key in selected],
                    "gpmag": [rows[key][1] for key in selected],
                    "rpmag": [rows[key][2] for key in selected],
                    "ipmag": [rows[key][3] for key in selected],
                    "zpmag": [rows[key][4] for key in selected],
                    "e_upmag": [0.01] * len(selected),
                    "e_gpmag": [0.02] * len(selected),
                    "e_rpmag": [0.03] * len(selected),
                    "e_ipmag": [0.04] * len(selected),
                    "e_zpmag": [0.05] * len(selected),
                }
            )
        ]


def test_enrich_sdss_psf_photometry_batches_objid_lookup():
    table = Table(
        {
            "source_id": ["a", "b", "c"],
            "objID": [123, 456, 789],
            "umag": [20.0, 21.0, 22.0],
        }
    )
    FakeVizier.calls = []

    enriched = enrich_sdss_psf_photometry(table, chunk_size=2, cache=False, vizier_cls=FakeVizier)

    assert FakeVizier.calls[0][1] == "123,456"
    assert FakeVizier.calls[1][1] == "789"
    assert enriched["upmag"].tolist()[:2] == [18.0, 19.0]
    assert enriched["gpmag"].tolist()[:2] == [17.0, 18.0]
    assert enriched["e_zpmag"].tolist()[:2] == [0.05, 0.05]
    assert np.isnan(enriched["upmag"][2])
    assert enriched["umag"].tolist() == [20.0, 21.0, 22.0]
