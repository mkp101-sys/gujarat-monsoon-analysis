"""Compare a fresh pipeline run with the reference run (results/reference/colab_run_analysis_outputs.zip).

Usage:
    python src/verify_outputs.py <new_output_folder> [--reference results/reference/colab_run_analysis_outputs.zip]

Checks that:
  1. every reference file exists in the new run (PNG and CSV, same relative paths);
  2. every CSV has the same columns and row labels;
  3. every numeric value agrees within a tolerance (default 0.05, which covers rounding
     differences between platforms);
  4. PNG files exist and are non-empty.
Exit code 0 = pass, 1 = fail.
"""
import argparse, io, os, sys, zipfile
import numpy as np
import pandas as pd

DEFAULT_REF = os.path.join(os.path.dirname(__file__), "..", "results", "reference", "colab_run_analysis_outputs.zip")

def read_reference(path):
    with zipfile.ZipFile(path) as z:
        return {n: z.read(n) for n in z.namelist() if not n.endswith("/")}

def compare_csv(name, ref_bytes, new_bytes, tol):
    problems = []
    ref = pd.read_csv(io.BytesIO(ref_bytes), index_col=0)
    new = pd.read_csv(io.BytesIO(new_bytes), index_col=0)
    if list(ref.columns) != list(new.columns):
        problems.append("columns differ")
    if list(map(str, ref.index)) != list(map(str, new.index)):
        problems.append("row labels differ")
        return problems
    for col in ref.columns:
        a = pd.to_numeric(ref[col], errors="coerce")
        b = pd.to_numeric(new[col], errors="coerce")
        if a.notna().any():
            diff = np.nanmax(np.abs(a.values - b.values))
            if not np.isfinite(diff) or diff > tol:
                problems.append(f"column '{col}' max difference {diff:.4f} > {tol}")
        elif not ref[col].astype(str).equals(new[col].astype(str)):
            problems.append(f"column '{col}' text differs")
    return problems

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("new_dir", help="folder produced by a new run (analysis_outputs/)")
    ap.add_argument("--reference", default=DEFAULT_REF)
    ap.add_argument("--tol", type=float, default=0.05)
    args = ap.parse_args()

    reference = read_reference(args.reference)
    failures = []
    checked = {"csv": 0, "png": 0, "other": 0}

    for name, data in reference.items():
        new_path = os.path.join(args.new_dir, name)
        if not os.path.exists(new_path):
            failures.append(f"missing: {name}")
            continue
        if name.endswith(".csv"):
            checked["csv"] += 1
            with open(new_path, "rb") as f:
                probs = compare_csv(name, data, f.read(), args.tol)
            failures += [f"{name}: {p}" for p in probs]
        elif name.endswith(".png"):
            checked["png"] += 1
            if os.path.getsize(new_path) == 0:
                failures.append(f"empty image: {name}")
        else:
            checked["other"] += 1

    print(f"Checked {checked['csv']} CSV, {checked['png']} PNG, {checked['other']} other files against the reference.")
    if failures:
        print(f"FAIL: {len(failures)} problem(s)")
        for f in failures[:50]:
            print("  -", f)
        sys.exit(1)
    print("PASS: the new run matches the reference within tolerance.")

if __name__ == "__main__":
    main()
