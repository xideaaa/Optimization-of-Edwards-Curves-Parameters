# Optimization-of-Edwards-Curves-Parameters
Optimization of Edwards Curves' Parameters over Finite Fields for Post-Quantum Isogeny-based Cryptography

------------------------------------------------------------
 Curve generation
------------------------------------------------------------

make gen N=1000 P=<hex> A=<hex> [COMPLETE=1] [OUT=path]

Generates N candidate twisted Edwards curves over the field F_p
with parameters a and d.

- If COMPLETE=1, only complete curves are generated
- Output is written to OUT in CSV format

------------------------------------------------------------
 Curve counting (SageMath)
------------------------------------------------------------

make count-file IN=<candidates.csv> OUT=<counted.csv>
[SMALL_BOUND=2 PRP=1 SKIP_TWIST=0 ALGO=sea]

Computes group orders using SageMath.

Environment variables passed to the Sage script:

  SMALL_BOUND  bound for stripping odd factors from the group order
  PRP          use pseudoprime test (1, default) or deterministic test (0)
  SKIP_TWIST   skip twist order computation if set to 1
  ALGO         point-counting algorithm (e.g. sea)

------------------------------------------------------------
 Curve selection
------------------------------------------------------------

make select-file IN=<counted.csv> OUT=<accepted.csv>

Selects curves with the required group structure.

------------------------------------------------------------
 Group operation benchmarks
------------------------------------------------------------

make bench-curves COST=<0/1> IN=<accepted.csv> OUT=<bench.csv>

Benchmarks group operations and scalar multiplication.
- If COST=1, arithmetic operation counts are recorded

------------------------------------------------------------
 Isogeny benchmarks
------------------------------------------------------------

make bench-iso COST=<0/1> IN=<accepted.csv> OUT=<bench_iso.csv>

Benchmarks isogeny-related operations.

After running bench-iso, execute the following scripts in order:
  python scripts/merge_one.py
  python scripts/merge_all.py
  python scripts/heuristic.py

These scripts merge benchmark results and compute heuristic scores
used to rank curve parameters.

------------------------------------------------------------
 Notes
------------------------------------------------------------
- SageMath is required for the counting stage
- All intermediate data is exchanged via CSV files
- The pipeline is intended to be run sequentially
