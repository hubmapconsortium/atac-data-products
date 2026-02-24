#!/usr/bin/env python3
"""
Quick check for genome coordinates in each modality of an .h5mu file.

Usage:
  python check_h5mu_bins.py --h5mu /data/.../HT_raw.h5mu
"""

from __future__ import annotations

import argparse
import re
import sys
from typing import Optional


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5mu", required=True, help="Path to input .h5mu")
    return p.parse_args()


def clean_chr(x: str) -> str:
    x = str(x).strip()
    x = re.sub(r"^chrchr", "chr", x)
    if x.startswith("chr"):
        return x
    if x.isdigit():
        return f"chr{x}"
    m = re.match(r"(\d+)", x)
    if m:
        return f"chr{m.group(1)}"
    return "chrUn"


def parse_var_name(name: str) -> Optional[tuple[str, str, Optional[str]]]:
    """Return (chrom, start, end_or_none) if the name looks like a bin."""
    name = str(name).strip()
    m = re.match(r"^([^:\s]+)[:\s]+([0-9eE.+-]+)[-:]+([0-9eE.+-]+)", name)
    if m:
        return clean_chr(m.group(1)), m.group(2), m.group(3)
    tokens = name.replace(",", " ").split()
    if len(tokens) >= 3:
        chrom_tok, start_tok, end_tok = tokens[-3], tokens[-2], tokens[-1]
        return clean_chr(chrom_tok), start_tok, end_tok
    if len(tokens) == 2:
        chrom_tok, start_tok = tokens
        return clean_chr(chrom_tok), start_tok, None
    return None


def main():
    args = parse_args()

    try:
        import muon as mu  # type: ignore
    except Exception as e:  # pragma: no cover
        sys.exit(f"muon is required to read .h5mu files; import failed: {e}")

    print(f"Loading {args.h5mu} ...")
    mdata = mu.read_h5mu(args.h5mu)
    print(f"Found modalities: {list(mdata.mod.keys())}")

    any_bins = False
    for key, ad in mdata.mod.items():
        var = ad.var
        print(f"\nChecking modality '{key}' with {ad.n_vars} vars ...")

        if "interval" in var.columns:
            print("  ✓ interval column present (chrom/start/end).")
            any_bins = True
            continue
        if "chrom" in var.columns and "start" in var.columns:
            print("  ✓ chrom/start columns present (end inferred if missing).")
            any_bins = True
            continue

        sample = [str(x) for x in ad.var_names[: min(50, ad.n_vars)]]
        parsed = [parse_var_name(s) for s in sample]
        if any(parsed):
            print("  ✓ var_names look like genomic bins (examples):")
            for s, p in zip(sample[:5], parsed[:5]):
                if p:
                    chrom, start, end = p
                    end_part = end if end is not None else "(none)"
                    print(f"    {s!r} -> {chrom}:{start}-{end_part}")
                    break
            any_bins = True
        else:
            print("  ✗ no obvious bin coordinates (var_names look like genes?)")

    if not any_bins:
        sys.exit("No modality with parseable genome coordinates found.")


if __name__ == "__main__":
    main()
