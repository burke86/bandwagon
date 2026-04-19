# bandwagon

Bulk cross-match helpers for collecting archival multi-band photometry, redshifts, and spectra from VizieR/CDS
catalogs.

## Installation

Install from GitHub:

```bash
python -m pip install "bandwagon @ git+https://github.com/burke86/bandwagon.git"
```

For local development:

```bash
git clone git@github.com:burke86/bandwagon.git
cd bandwagon
python -m pip install -e ".[dev]"
```

You can also install the pinned runtime list directly:

```bash
python -m pip install -r requirements.txt
```

Run the tests with:

```bash
pytest -q
```

## Usage

```python
from astropy.coordinates import SkyCoord
import astropy.units as u
from bandwagon import (
    DEFAULT_CATALOGS,
    matches_to_photometry,
    query_archival_spectra,
    query_simbad_redshifts,
    select_best_redshift,
    xmatch_catalogs,
)

coords = SkyCoord(
    ra=[10.0, 11.0] * u.deg,
    dec=[-2.0, -3.0] * u.deg,
    frame="icrs",
)

matches = xmatch_catalogs(coords, source_id=["src-a", "src-b"])
photometry = matches_to_photometry(matches)

redshift_candidates = query_simbad_redshifts(coords, source_id=["src-a", "src-b"])
redshifts = select_best_redshift(redshift_candidates)

spectra = query_archival_spectra(
    coords,
    source_id=["src-a", "src-b"],
    providers=("desi", "sdss"),
)
```

SIMBAD redshift queries return a separate long-form table with `source_id`,
`object_name`, `redshift`, `redshift_err`, `quality`, `reference`, and
`match_distance_arcsec`. They are not folded into the photometry table, so
conflicting redshift candidates can be inspected before selecting one.

Spectra queries are opt-in and return metadata/index rows; they do not download
spectra by default. Supported providers are `desi`, `sdss`, `lamost`, `6dfgs`,
and `mast`. DESI uses SPARCL and requires the optional spectra extra:

```bash
python -m pip install -e ".[spectra]"
```

LAMOST and 6dFGS are metadata-first VizieR XMatch providers. SDSS and MAST use
`astroquery`. Bandwagon does not use `pyvo`.

The default catalog set is:

| Output key | VizieR table | Bands | Match radius |
| --- | --- | --- | ---: |
| `galex_ais` | `II/335/galex_ais` | `FUV`, `NUV` | `3.0"` |
| `sdss_dr16` | `V/154/sdss16` | `u`, `g`, `r`, `i`, `z` | `1.0"` |
| `2mass` | `II/246/out` | `J`, `H`, `Ks` | `2.0"` |
| `allwise` | `II/328/allwise` | `W1`, `W2`, `W3`, `W4` | `3.0"` |

Optional catalog aliases are also available:

| Alias | VizieR table | Bands | Default radius |
| --- | --- | --- | ---: |
| `akari_irc` | `II/297/irc` | `S9W`, `L18W` | `6.0"` |
| `akari_fis` | `II/298/fis` | `N60`, `WIDE-S`, `WIDE-L`, `N160` | `20.0"` |
| `iras_psc` | `II/125/main` | `F12`, `F25`, `F60`, `F100` | `30.0"` |

```python
matches = xmatch_catalogs(
    coords,
    catalogs={
        **DEFAULT_CATALOGS,
        "2mass": "2mass",
        "akari_irc": "akari_irc",
        "akari_fis": "akari_fis",
        "iras_psc": "iras_psc",
    },
    source_id=["src-a", "src-b"],
)
photometry = matches_to_photometry(matches, min_quality=2)
```

`2mass` is converted from Vega magnitudes. AKARI and IRAS publish flux
densities, so Bandwagon converts Jy to mJy directly. IRAS uncertainty columns
are percent flux uncertainties; AKARI uncertainty columns are Jy.

DESI/Legacy Survey aliases are available for VizieR's DR8 north/south
photometric redshift tables, but SDSS is the default optical source because it
provides the full `ugriz` set through CDS XMatch.
