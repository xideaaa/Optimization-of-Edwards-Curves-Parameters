#!/usr/bin/env python3

import pandas as pd
from pathlib import Path


DATA_DIR = Path("data")

REQUIRED = [
    "accepted.csv",
    "bench_time.csv",
    "bench_iso.csv"
]

def merge_one_dir(base: Path):
    acc_path  = base / "accepted.csv"
    time_path = base / "bench_time.csv"
    iso_path  = base / "bench_iso.csv"

    if not (acc_path.exists() and time_path.exists() and iso_path.exists()):
        return False

    print(f"[+] merging {base}")

    df_acc  = pd.read_csv(acc_path)
    df_time = pd.read_csv(time_path)
    df_iso  = pd.read_csv(iso_path)

    df_acc = df_acc[
        [
            "id",
            "order_hex",
            "h2",
            "L_hex",
            "bits_L",
            "twist_h2",
            "twist_L_hex",
            "twist_bits_L"
        ]
    ].copy()

    df_time = df_time[
        [
            "id","p_hex","a_hex","d_hex",
            "algo","win",
            "ns_avg","M_avg","S_avg","A_avg"
        ]
    ].copy()

    df_iso = df_iso[
        [
            "id","p_hex","a_hex","d_hex",
            "stage","k",
            "ns_avg","M_avg","S_avg","A_avg"
        ]
    ].copy()

    iso_pivot = (
        df_iso
        .pivot_table(
            index=["id","p_hex","a_hex","d_hex"],
            columns=["stage","k"],
            values=["ns_avg","M_avg","S_avg","A_avg"]
        )
    )

    iso_pivot.columns = [
        f"iso_{metric}_{stage}_{k}"
        for metric,stage,k in iso_pivot.columns
    ]

    iso_pivot = iso_pivot.reset_index()

    time_pivot = (
        df_time
        .pivot_table(
            index=["id","p_hex","a_hex","d_hex"],
            columns=["algo","win"],
            values=["ns_avg","M_avg","S_avg","A_avg"]
        )
    )

    time_pivot.columns = [
        f"smul_{metric}_{algo}_w{win}"
        for metric,algo,win in time_pivot.columns
    ]

    time_pivot = time_pivot.reset_index()

    merged = pd.merge(
        time_pivot,
        iso_pivot,
        on=["id","p_hex","a_hex","d_hex"],
        how="inner"
    )

    merged = pd.merge(
        merged,
        df_acc,
        on="id",
        how="left"
    )

    out = base / "bench_merged.csv"
    merged.to_csv(out, index=False)

    print(f"    -> wrote {out}")
    return True

def main():
    count = 0
    for sub in sorted(DATA_DIR.iterdir()):
        if sub.is_dir():
            if merge_one_dir(sub):
                count += 1

if __name__ == "__main__":
    main()
