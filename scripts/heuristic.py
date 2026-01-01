#!/usr/bin/env python3
import math
import random
import pandas as pd
from copy import deepcopy

def normalize(series):
    return (series - series.min()) / (series.max() - series.min() + 1e-12)

def clamp(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))

def prepare_dataframe(df):
    """
    Uses already-aggregated columns produced by merge scripts.
    """

    # curve arithmetic cost
    df["curve_cost_raw"] = (
        df["smul_ns_avg_ladder_w0"] +
        df["smul_ns_avg_add_w0"] +
        df["smul_ns_avg_double_w0"]
    )

    # isogeny cost
    df["iso_cost_raw"] = (
        df["iso_ns_avg_kernel_1"] +
        df["iso_ns_avg_eval_1"] +
        df["iso_ns_avg_chain_10"] +
        df["iso_ns_avg_chain_20"] +
        df["iso_ns_avg_chain_50"]
    )

    # omega penalty
    df["omega_cost_raw"] = (
        -df["h2"].apply(lambda x: math.log2(x)) -
        df["twist_h2"].apply(lambda x: math.log2(x)) -
        df["bits_L"]
    )

    df["curve_cost"] = normalize(df["curve_cost_raw"])
    df["iso_cost"]   = normalize(df["iso_cost_raw"])
    df["omega_cost"] = normalize(df["omega_cost_raw"])

    return df

def objective(df, w):
    if w["alpha"] < 0.05 or w["beta"] < 0.05 or w["gamma"] < 0.05:
        return 1e9

    score = (
        w["alpha"] * df["curve_cost"] +
        w["beta"]  * df["iso_cost"] +
        w["gamma"] * df["omega_cost"]
    )
    return score.min()

def random_weights():
    a = random.random()
    b = random.random()
    c = random.random()
    s = a + b + c
    return {"alpha": a/s, "beta": b/s, "gamma": c/s}

def neighbor(w, scale=0.1):
    w2 = deepcopy(w)
    k = random.choice(["alpha", "beta", "gamma"])
    w2[k] += random.uniform(-scale, scale)

    for x in ["alpha","beta","gamma"]:
        w2[x] = clamp(w2[x])

    s = w2["alpha"] + w2["beta"] + w2["gamma"]
    if s == 0:
        return random_weights()

    for x in ["alpha","beta","gamma"]:
        w2[x] /= s

    return w2

# random search
def random_search(df, iters=5000):
    best_w = None
    best_s = float("inf")
    for _ in range(iters):
        w = random_weights()
        s = objective(df, w)
        if s < best_s:
            best_s = s
            best_w = w
    return best_w, best_s

# simulated annealing
def simulated_annealing(df, w0, T0=1.0, Tmin=1e-4, alpha=0.995, steps=10000):
    w = deepcopy(w0)
    s = objective(df, w)

    best_w = deepcopy(w)
    best_s = s
    T = T0

    for _ in range(steps):
        w2 = neighbor(w)
        s2 = objective(df, w2)

        if s2 < s or random.random() < math.exp((s - s2)/T):
            w, s = w2, s2
            if s < best_s:
                best_w, best_s = deepcopy(w), s

        T *= alpha
        if T < Tmin:
            break

    return best_w, best_s

def main():
    in_csv  = "data/heuristic_input.csv"
    out_csv = "data/heuristic_scored.csv"
    weights_csv = "data/heuristic_weights_runs.csv"
    top_csv = "data/heuristic_top10.csv"

    df = pd.read_csv(in_csv)
    df = prepare_dataframe(df)

    records = []
    N_RUNS = 1000

    for i in range(N_RUNS):
        random.seed(i)
        w0, _ = random_search(df)
        w, score = simulated_annealing(df, w0)

        records.append({
            "run": i,
            "alpha": w["alpha"],
            "beta":  w["beta"],
            "gamma": w["gamma"],
            "best_score": score
        })

    weights_df = pd.DataFrame(records)
    weights_df.to_csv(weights_csv, index=False)

    w_med = {
        "alpha": weights_df["alpha"].median(),
        "beta":  weights_df["beta"].median(),
        "gamma": weights_df["gamma"].median()
    }

    s = sum(w_med.values())
    for k in w_med:
        w_med[k] /= s

    df["H_final"] = (
        w_med["alpha"] * df["curve_cost"] +
        w_med["beta"]  * df["iso_cost"] +
        w_med["gamma"] * df["omega_cost"]
    )

    df.sort_values("H_final", inplace=True)
    df.to_csv(out_csv, index=False)

    top10_best  = df.head(10).copy()
    top10_worst = df.tail(10).copy()

    top10_best["rank"] = range(1, 11)
    top10_worst["rank"] = range(len(df)-9, len(df)+1)

    extremes = pd.concat([top10_best, top10_worst])
    extremes.to_csv(top_csv, index=False)

if __name__ == "__main__":
    main()
