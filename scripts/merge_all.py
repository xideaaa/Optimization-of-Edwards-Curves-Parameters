#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

DATA_DIR = Path("data")
OUT_FILE = DATA_DIR / "heuristic_input.csv"

def main():
    rows = []
    used_dirs = 0

    for sub in sorted(DATA_DIR.iterdir()):
        if not sub.is_dir():
            continue

        f = sub / "bench_merged.csv"
        if not f.exists():
            continue

        df = df.copy()
        df["curve_dir"] = sub.name

        rows.append(df)
        used_dirs += 1

    merged = pd.concat(rows, ignore_index=True)
    merged.to_csv(OUT_FILE, index=False)

if __name__ == "__main__":
    main()
