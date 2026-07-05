#!/usr/bin/env python
"""Build chunked Stripe 82 parent catalogs from Legacy Survey sweep files.

This script expects local Legacy Survey DR10.1 sweep FITS files, filters them to
the Stripe 82 footprint, adds a stable source_id, and writes HDF5 chunks.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import urllib.request
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack


DEFAULT_COLUMNS = (
    "RELEASE",
    "BRICKID",
    "BRICKNAME",
    "OBJID",
    "TYPE",
    "RA",
    "DEC",
    "EBV",
    "MASKBITS",
    "FLUX_G",
    "FLUX_R",
    "FLUX_Z",
    "FLUX_W1",
    "FLUX_W2",
    "FLUX_W3",
    "FLUX_W4",
    "FLUX_IVAR_G",
    "FLUX_IVAR_R",
    "FLUX_IVAR_Z",
    "FLUX_IVAR_W1",
    "FLUX_IVAR_W2",
    "FLUX_IVAR_W3",
    "FLUX_IVAR_W4",
    "MW_TRANSMISSION_G",
    "MW_TRANSMISSION_R",
    "MW_TRANSMISSION_Z",
    "MW_TRANSMISSION_W1",
    "MW_TRANSMISSION_W2",
    "MW_TRANSMISSION_W3",
    "MW_TRANSMISSION_W4",
    "NOBS_G",
    "NOBS_R",
    "NOBS_Z",
    "NOBS_W1",
    "NOBS_W2",
    "FRACFLUX_G",
    "FRACFLUX_R",
    "FRACFLUX_Z",
    "FRACFLUX_W1",
    "FRACFLUX_W2",
    "FRACMASKED_G",
    "FRACMASKED_R",
    "FRACMASKED_Z",
    "FRACIN_G",
    "FRACIN_R",
    "FRACIN_Z",
    "ANYMASK_G",
    "ANYMASK_R",
    "ANYMASK_Z",
    "ALLMASK_G",
    "ALLMASK_R",
    "ALLMASK_Z",
    "WISEMASK_W1",
    "WISEMASK_W2",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter Legacy Survey sweep FITS files to Stripe 82 and write HDF5 chunks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "sweep_dir",
        type=Path,
        help="Directory containing Legacy Survey sweep-*.fits files.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download missing DR10.1 sweep files overlapping the requested Stripe 82 footprint.",
    )
    parser.add_argument(
        "--base-url",
        default="https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/sweep/10.1",
        help="Base URL for Legacy Survey sweep files when --download is set.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("stripe82_legacy_parent_chunks"),
        help="Directory for HDF5 chunk files and manifest.csv.",
    )
    parser.add_argument(
        "--glob",
        default="sweep-*.fits",
        help="File glob to select sweep files within sweep_dir.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100_000,
        help="Maximum number of Stripe 82 objects per output HDF5 chunk.",
    )
    parser.add_argument(
        "--hdf5-path",
        default="parent",
        help="Internal HDF5 dataset path used by astropy Table.write.",
    )
    parser.add_argument(
        "--prefix",
        default="stripe82_legacy_parent",
        help="Output filename prefix for chunk files.",
    )
    parser.add_argument(
        "--dec-min",
        type=float,
        default=-1.266,
        help="Minimum Stripe 82 declination in degrees.",
    )
    parser.add_argument(
        "--dec-max",
        type=float,
        default=1.266,
        help="Maximum Stripe 82 declination in degrees.",
    )
    parser.add_argument(
        "--ra-west-min",
        type=float,
        default=300.0,
        help="Western Stripe 82 RA lower bound in degrees.",
    )
    parser.add_argument(
        "--ra-east-max",
        type=float,
        default=60.0,
        help="Eastern Stripe 82 RA upper bound in degrees.",
    )
    parser.add_argument(
        "--columns",
        choices=("sed", "all"),
        default="sed",
        help="Write a compact SED-oriented column subset or every sweep column.",
    )
    parser.add_argument(
        "--no-source-id",
        action="store_true",
        help="Do not add source_id=RELEASE-BRICKID-OBJID.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output chunk files and manifest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive.")

    if args.download:
        download_sweeps(
            args.sweep_dir,
            base_url=args.base_url,
            dec_min=args.dec_min,
            dec_max=args.dec_max,
            ra_west_min=args.ra_west_min,
            ra_east_max=args.ra_east_max,
        )

    sweep_paths = sorted(args.sweep_dir.glob(args.glob))
    if not sweep_paths:
        raise FileNotFoundError(f"No files matching {args.glob!r} under {args.sweep_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "manifest.csv"
    if manifest_path.exists() and not args.overwrite:
        raise FileExistsError(f"{manifest_path} exists; pass --overwrite to replace it.")

    pending: list[Table] = []
    pending_rows = 0
    chunk_index = 0
    total_rows = 0
    manifest_rows = []

    for sweep_path in sweep_paths:
        table = Table.read(sweep_path, memmap=True)
        filtered = filter_stripe82(
            table,
            dec_min=args.dec_min,
            dec_max=args.dec_max,
            ra_west_min=args.ra_west_min,
            ra_east_max=args.ra_east_max,
        )
        if len(filtered) == 0:
            continue

        if args.columns == "sed":
            filtered = keep_existing_columns(filtered, DEFAULT_COLUMNS)
        if not args.no_source_id:
            add_source_id(filtered)

        start = 0
        while start < len(filtered):
            take = min(args.chunk_size - pending_rows, len(filtered) - start)
            pending.append(filtered[start : start + take])
            pending_rows += take
            start += take

            if pending_rows == args.chunk_size:
                row = write_chunk(
                    pending,
                    output_dir=args.output_dir,
                    prefix=args.prefix,
                    chunk_index=chunk_index,
                    hdf5_path=args.hdf5_path,
                    overwrite=args.overwrite,
                )
                manifest_rows.append(row)
                total_rows += row["n_rows"]
                chunk_index += 1
                pending = []
                pending_rows = 0

    if pending_rows:
        row = write_chunk(
            pending,
            output_dir=args.output_dir,
            prefix=args.prefix,
            chunk_index=chunk_index,
            hdf5_path=args.hdf5_path,
            overwrite=args.overwrite,
        )
        manifest_rows.append(row)
        total_rows += row["n_rows"]

    write_manifest(manifest_path, manifest_rows, overwrite=args.overwrite)
    print(f"Wrote {len(manifest_rows)} chunks with {total_rows} Stripe 82 objects.")
    print(f"Manifest: {manifest_path}")


def download_sweeps(
    sweep_dir: Path,
    *,
    base_url: str,
    dec_min: float,
    dec_max: float,
    ra_west_min: float,
    ra_east_max: float,
) -> None:
    sweep_dir.mkdir(parents=True, exist_ok=True)
    names = sweep_names_for_footprint(
        dec_min=dec_min,
        dec_max=dec_max,
        ra_west_min=ra_west_min,
        ra_east_max=ra_east_max,
    )
    base_url = base_url.rstrip("/")
    for i, name in enumerate(names, start=1):
        target = sweep_dir / name
        if target.exists():
            print(f"[{i}/{len(names)}] exists {target}")
            continue
        url = f"{base_url}/{name}"
        tmp = target.with_suffix(target.suffix + ".part")
        print(f"[{i}/{len(names)}] downloading {url}")
        try:
            with urllib.request.urlopen(url) as response, tmp.open("wb") as handle:
                shutil.copyfileobj(response, handle)
            tmp.replace(target)
        except Exception:
            tmp.unlink(missing_ok=True)
            raise


def sweep_names_for_footprint(
    *,
    dec_min: float,
    dec_max: float,
    ra_west_min: float,
    ra_east_max: float,
) -> list[str]:
    names = []
    for ra1 in range(0, 360, 5):
        ra2 = ra1 + 5
        if not ra_interval_overlaps_wrapped_footprint(
            ra1,
            ra2,
            ra_west_min=ra_west_min,
            ra_east_max=ra_east_max,
        ):
            continue
        for dec1 in range(-90, 40, 5):
            dec2 = dec1 + 5
            if dec2 < dec_min or dec1 > dec_max:
                continue
            names.append(f"sweep-{sky_token(ra1, dec1)}-{sky_token(ra2, dec2)}.fits")
    return names


def ra_interval_overlaps_wrapped_footprint(
    ra1: int,
    ra2: int,
    *,
    ra_west_min: float,
    ra_east_max: float,
) -> bool:
    return (ra2 > ra_west_min and ra1 < 360.0) or (ra1 < ra_east_max and ra2 > 0.0)


def sky_token(ra: int, dec: int) -> str:
    sign = "p" if dec >= 0 else "m"
    return f"{ra:03d}{sign}{abs(dec):03d}"


def filter_stripe82(
    table: Table,
    *,
    dec_min: float,
    dec_max: float,
    ra_west_min: float,
    ra_east_max: float,
) -> Table:
    ra_col = column_name(table, "RA")
    dec_col = column_name(table, "DEC")
    ra = np.asarray(table[ra_col], dtype=float)
    dec = np.asarray(table[dec_col], dtype=float)
    mask = (dec >= dec_min) & (dec <= dec_max) & ((ra >= ra_west_min) | (ra <= ra_east_max))
    return table[mask]


def keep_existing_columns(table: Table, requested: tuple[str, ...]) -> Table:
    lookup = {name.upper(): name for name in table.colnames}
    keep = [lookup[name] for name in requested if name in lookup]
    return table[keep]


def add_source_id(table: Table) -> None:
    release_col = column_name(table, "RELEASE")
    brickid_col = column_name(table, "BRICKID")
    objid_col = column_name(table, "OBJID")
    table["source_id"] = np.array(
        [
            f"{release}-{brickid}-{objid}"
            for release, brickid, objid in zip(
                table[release_col],
                table[brickid_col],
                table[objid_col],
                strict=True,
            )
        ],
        dtype="U64",
    )


def column_name(table: Table, target: str) -> str:
    target_upper = target.upper()
    for name in table.colnames:
        if name.upper() == target_upper:
            return name
    raise KeyError(f"Required column {target!r} not found.")


def write_chunk(
    tables: list[Table],
    *,
    output_dir: Path,
    prefix: str,
    chunk_index: int,
    hdf5_path: str,
    overwrite: bool,
) -> dict[str, object]:
    chunk = vstack(tables, metadata_conflicts="silent") if len(tables) > 1 else tables[0]
    output_path = output_dir / f"{prefix}_{chunk_index:05d}.h5"
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"{output_path} exists; pass --overwrite to replace it.")

    chunk.write(output_path, path=hdf5_path, format="hdf5", serialize_meta=True, overwrite=overwrite)

    ra_col = column_name(chunk, "RA")
    dec_col = column_name(chunk, "DEC")
    return {
        "chunk": chunk_index,
        "filename": output_path.name,
        "n_rows": len(chunk),
        "ra_min": float(np.nanmin(chunk[ra_col])),
        "ra_max": float(np.nanmax(chunk[ra_col])),
        "dec_min": float(np.nanmin(chunk[dec_col])),
        "dec_max": float(np.nanmax(chunk[dec_col])),
    }


def write_manifest(path: Path, rows: list[dict[str, object]], *, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass --overwrite to replace it.")
    fieldnames = ["chunk", "filename", "n_rows", "ra_min", "ra_max", "dec_min", "dec_max"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
