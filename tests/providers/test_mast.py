import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from bandwagon.providers.mast import query_mast_spectra


class FakeObservations:
    @staticmethod
    def query_region(coord, *, radius):
        return Table(
            {
                "obs_id": ["obs1", "obs2"],
                "dataproduct_type": ["spectrum", "image"],
            }
        )

    @staticmethod
    def get_product_list(obs):
        assert len(obs) == 1
        return Table(
            {
                "obs_id": ["obs1"],
                "obs_collection": ["HST"],
                "instrument_name": ["COS"],
                "target_name": ["Fairall 9"],
                "calib_level": [2],
                "dataURI": ["mast:product/file.fits"],
                "productSubGroupDescription": ["X1D"],
                "productFilename": ["file.fits"],
                "s_ra": [20.0],
                "s_dec": [-10.0],
            }
        )

    @staticmethod
    def filter_products(products, *, productType):
        assert productType == ["SCIENCE"]
        return products


def test_query_mast_spectra_keeps_spectrum_products():
    coords = SkyCoord([20.0] * u.deg, [-10.0] * u.deg)

    table = query_mast_spectra(coords, source_id=["src"], client=FakeObservations)

    assert len(table) == 1
    assert table["provider"][0] == "mast"
    assert table["survey"][0] == "HST"
    assert table["instrument"][0] == "COS"
    assert table["access_url"][0] == "mast:product/file.fits"
    assert table["mast_product_filename"][0] == "file.fits"
