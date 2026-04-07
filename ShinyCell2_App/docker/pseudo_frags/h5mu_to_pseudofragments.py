#!/usr/bin/env python3
"""
Convert an h5mu ATAC-by-bin modality into a bgzip-compressed, tabix-indexed
pseudo fragment file. Coordinates tolerate scientific notation (e.g. 1e+05).

Usage:
  python h5mu_to_pseudofragments_new.py --h5mu HT.h5mu --out pseudo_HT_fragments.tsv.gz
"""

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import muon as mu
import numpy as np
import scipy.sparse as sp


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Expand an ATAC-by-bin matrix into a pseudo fragment table."
    )
    parser.add_argument("--h5mu", required=True, help="Input .h5mu file")
    parser.add_argument(
        "--out",
        required=True,
        help="Output path (compressed with bgzip; .gz/.bgz added if missing)",
    )
    parser.add_argument(
        "--modality",
        default=None,
        help="Modality key; defaults to 'atac_cell_by_bin' or the first modality.",
    )
    return parser.parse_args()


def _require_cmd(name: str):
    if shutil.which(name) is None:
        raise RuntimeError(
            f"Required command '{name}' not found. Install htslib/samtools."
        )


def _bgzip_and_tabix(src: Path, dest: Path):
    """Compress src -> dest with bgzip and create a tabix index."""
    _require_cmd("bgzip")
    _require_cmd("tabix")

    print(f"Compressing {src} -> {dest} with bgzip ...")
    with dest.open("wb") as out_f:
        subprocess.run(["bgzip", "-c", str(src)], stdout=out_f, check=True)

    print(f"Indexing {dest} with tabix ...")
    subprocess.run(
        ["tabix", "-s", "1", "-b", "2", "-e", "3", "-0", "-f", str(dest)],
        check=True,
    )

    src.unlink(missing_ok=True)


# ------------------ chromosome + coordinate cleaners ------------------

def clean_chr(x: str) -> str:
    """Normalize chromosome labels robustly."""
    x = str(x).strip()

    # remove double-prefix
    x = re.sub(r"^chrchr", "chr", x)

    if x.startswith("chr"):
        return x

    # pure numeric → chr##
    if x.isdigit():
        return f"chr{x}"

    # capture leading number
    m = re.match(r"(\d+)", x)
    if m:
        return f"chr{m.group(1)}"

    return "chrUn"


def extract_int(x: str) -> int:
    """Extract safe integer start coordinate."""
    try:
        return max(0, int(float(x)))
    except Exception:
        m = re.search(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", str(x))
        if m:
            return max(0, int(float(m.group(0))))
        return 0


def _to_int_array(values):
    return np.array([extract_int(v) for v in values], dtype=np.int64)


def _parse_var_name(name: str):
    """Best-effort parsing of a var name into (chrom, start, end_or_None)."""
    name = str(name).strip()

    # common pattern: chr1:0-500 or 1:0-500 or chr1 0 500
    m = re.match(r"^([^:\s]+)[:\s]+([0-9eE.+-]+)[-:]+([0-9eE.+-]+)", name)
    if m:
        return clean_chr(m.group(1)), m.group(2), m.group(3)

    # whitespace separated; take the last three tokens to be robust to prefixes
    tokens = name.replace(",", " ").split()
    if len(tokens) >= 3:
        chrom_tok, start_tok, end_tok = tokens[-3], tokens[-2], tokens[-1]
        return clean_chr(chrom_tok), start_tok, end_tok

    if len(tokens) == 2:  # chrom + start only; synthesize end later
        chrom_tok, start_tok = tokens
        return clean_chr(chrom_tok), start_tok, None

    return None


# ----------------------------------------------------------------------
# Modality and bin loading
# ----------------------------------------------------------------------

def pick_modality(mdata, override: str | None):
    if override:
        if override in mdata.mod:
            return mdata.mod[override]
        raise ValueError(f"Modality '{override}' not found.")
    if "atac_cell_by_bin" in mdata.mod:
        return mdata.mod["atac_cell_by_bin"]

    # try every modality until bin parsing succeeds
    errors: dict[str, str] = {}
    for key, ad in mdata.mod.items():
        try:
            load_correct_bins(ad)
            print(f"Auto-selected modality '{key}' with parseable bin coordinates.")
            return ad
        except Exception as e:
            errors[key] = str(e)

    available = ", ".join(mdata.mod.keys())
    raise ValueError(
        "Could not find a modality with parseable bin coordinates. "
        "Pass --modality to pick the correct ATAC/bin modality. "
        f"Available modalities: {available}. "
        f"Last errors: {errors}"
    )


def load_correct_bins(adata):
    var = adata.var

    # Case 1: interval column exists
    if "interval" in var.columns:
        intervals = var["interval"].astype(str).to_numpy()
        parts = np.array([s.split() for s in intervals])
        chrom = np.array([clean_chr(x) for x in parts[:, 0]], dtype=object)
        start = _to_int_array(parts[:, 1])
        end = _to_int_array(parts[:, 2])
        return chrom, end-500, end

    # Case 2: explicit chrom/start/end columns
    if "chrom" in var.columns and "start" in var.columns:
        chrom = np.array([clean_chr(x) for x in var["chrom"].astype(str)], dtype=object)
        start = _to_int_array(var["start"].to_numpy())
        if "end" in var.columns:
            end = _to_int_array(var["end"].to_numpy())
        else:
            end = start + 500
        return chrom, end-500, end

    # Case 3: fallback parse var_names
    parsed = [_parse_var_name(str(name)) for name in adata.var_names]
    if any(p is None for p in parsed):
        bad_idx = [i for i, p in enumerate(parsed) if p is None][0]
        bad = adata.var_names[bad_idx]
        raise ValueError(
            "Cannot parse chrom/start/end from var_names. "
            "Example failing entry: "
            f"{bad!r}. "
            "Pass --modality to select the correct ATAC/bin modality."
        )

    chrom = np.array([p[0] for p in parsed], dtype=object)
    start = _to_int_array([p[1] for p in parsed])
    raw_end = [p[2] for p in parsed]
    if any(e is None for e in raw_end):
        end = start + 500
    else:
        end = _to_int_array(raw_end)

    return chrom, end-500, end


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    args = parse_args()

    h5mu_path = Path(args.h5mu)
    out_path = Path(args.out)

    # enforce .gz or .bgz
    if out_path.suffix not in {".gz", ".bgz"}:
        gz_path = out_path.with_suffix(out_path.suffix + ".gz")
        print(f"Output will be bgzip-compressed as {gz_path}")
    else:
        gz_path = out_path

    gz_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = gz_path.with_name(gz_path.name + ".tmp")

    print(f"Loading MuData from {h5mu_path} ...")
    mdata = mu.read_h5mu(h5mu_path)
    adata = pick_modality(mdata, args.modality)
    print(f"Selected modality with {adata.n_obs} cells and {adata.n_vars} bins.")

    # matrix
    mat = adata.layers.get("counts", adata.X)
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)

    coo = mat.tocoo()

    chrom, start, end = load_correct_bins(adata)
    cells = adata.obs_names.to_numpy()

    print("Building fragment rows in memory for sorting ...")
    unsorted_path = gz_path.with_name(gz_path.name + ".unsorted")

    print(f"Writing unsorted fragments → {unsorted_path}")
    row_count = 0
    with unsorted_path.open("w") as f:
        for r, c, v in zip(coo.row, coo.col, coo.data):
            v = int(v)
            if v > 0:
                f.write(f"{chrom[c]}\t{start[c]}\t{end[c]}\t{cells[r]}\t{v}\n")
                row_count += 1

    print("Sorting fragments by chrom/start/end for tabix using `sort` ...")
    _require_cmd("sort")
    with tmp_path.open("w") as sorted_f:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", "-k3,3n", str(unsorted_path)],
            stdout=sorted_f,
            check=True,
        )

    unsorted_path.unlink(missing_ok=True)
    print(f"Wrote and sorted {row_count} fragments → {tmp_path}")

    print("Compressing + indexing ...")
    _bgzip_and_tabix(tmp_path, gz_path)

    print(f"Done: {gz_path} (index: {gz_path}.tbi)")


if __name__ == "__main__":
    main()
