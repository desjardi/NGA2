#!/usr/bin/env python3
"""Compare amr_shocktube output against a reference solution."""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


COLUMN_ALIASES = {
    "x":   ["position", "x", "x_cc", "xpos", "var1"],
    "vf":  ["vf", "volume_fraction", "alpha", "alphaL", "alpha_l"],
    "rho": ["rho", "density", "mixrho", "rho_mix", "var2"],
    "p":   ["p", "pressure", "press", "var3"],
    "u":   ["velocity", "u", "ux", "vel", "var4"],
}


def find_col(df, key):
    cols_lower = {c.strip().lower(): c for c in df.columns}
    for alias in COLUMN_ALIASES[key]:
        if alias in cols_lower:
            return cols_lower[alias]
    return None


def load(path):
    df = pd.read_csv(path)
    out = {}
    for k in COLUMN_ALIASES:
        c = find_col(df, k)
        if c is not None:
            out[k] = df[c].to_numpy()
    if "x" not in out:
        sys.exit(f"[error] no position column in {path}; columns are {list(df.columns)}")
    order = np.argsort(out["x"])
    return {k: v[order] for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", default=None,
                    help="Sim CSV. Defaults to the latest data/data_*.csv next to this script.")
    ap.add_argument("--ref", default=str(Path(__file__).resolve().parent / "shocktube_ref.csv"))
    ap.add_argument("--old", default=str(Path(__file__).resolve().parent / "oldcode.csv"))
    ap.add_argument("--out", default="compare.png")
    args = ap.parse_args()

    if args.sim is None:
        here = Path(__file__).resolve().parent
        candidates = sorted((here / "data").glob("data_*.csv"))
        if not candidates:
            sys.exit(f"[error] no data/data_*.csv in {here}; pass --sim explicitly")
        args.sim = str(candidates[-1])

    sim = load(args.sim)
    ref = load(args.ref)
    old = load(args.old) if Path(args.old).exists() else None

    panels = [("vf",  "Volume fraction"),
              ("rho", "Density"),
              ("p",   "Pressure"),
              ("u",   "Velocity")]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, (key, title) in zip(axes.flat, panels):
        if key in ref:
            ax.plot(ref["x"], ref[key], "k-", lw=1.2, label="reference")
        if old is not None and key in old:
            ax.plot(old["x"], old[key], "C0s", ms=3, mfc="none", label="old code")
        if key in sim:
            ax.plot(sim["x"], sim[key], "C3o", ms=3, mfc="none", label="amr_shocktube")
        ax.set_title(title)
        ax.set_xlabel("x")
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=9)
    suptitle = f"sim: {Path(args.sim).name}    ref: {Path(args.ref).name}"
    if old is not None:
        suptitle += f"    old: {Path(args.old).name}"
    fig.suptitle(suptitle)
    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")
    plt.show()


if __name__ == "__main__":
    main()
