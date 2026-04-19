# bandwagon

Bulk cross-match helpers for collecting multi-band photometry from VizieR/CDS
catalogs.

```python
from astropy.coordinates import SkyCoord
import astropy.units as u
from bandwagon import DEFAULT_CATALOGS, matches_to_photometry, xmatch_catalogs

coords = SkyCoord(
    ra=[10.0, 11.0] * u.deg,
    dec=[-2.0, -3.0] * u.deg,
    frame="icrs",
)

matches = xmatch_catalogs(coords, source_id=["src-a", "src-b"])
photometry = matches_to_photometry(matches)
```

The default catalog set is:

| Output key | VizieR table | Bands |
| --- | --- | --- |
| `galex_ais` | `II/335/galex_ais` | `FUV`, `NUV` |
| `sdss_dr16` | `V/154/sdss16` | `u`, `g`, `r`, `i`, `z` |
| `allwise` | `II/328/allwise` | `W1`, `W2`, `W3`, `W4` |

Optional catalog aliases are also available:

| Alias | VizieR table | Bands |
| --- | --- | --- |
| `2mass` | `II/246/out` | `J`, `H`, `Ks` |
| `akari_irc` | `II/297/irc` | `S9W`, `L18W` |
| `akari_fis` | `II/298/fis` | `N60`, `WIDE-S`, `WIDE-L`, `N160` |
| `iras_psc` | `II/125/main` | `F12`, `F25`, `F60`, `F100` |

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
