#!/usr/bin/env sage
# -*- coding: utf-8 -*-
#
# Counting points on twisted Edwards curves:
#   - read a CSV file with curve candidates (p, a, d)
#   - construct elliptic curve E(Fp) isomorphic to twisted Edwards
#   - compute #E
#   - factor #E as 2^v2 * L
#
# Flags:
#   PRP=1 - use is_pseudoprime (fast, default)
#   PRP=0 - use is_prime(proof=True) (slower, deterministic)
#   SKIP_TWIST=1 - skip twist cardinality
#   SMALL_BOUND=100000 - remove odd prime factors of L up to this bound
#
# Input (candidates.csv)
# Header:
#   id,p_hex,a_hex,d_hex,a_is_square,d_is_square,complete_ok
#
# Output (counted.csv):
#   id,p_hex,a_hex,d_hex,order_hex,h2,L_hex,is_L_prime,bits_L,
#   twist_order_hex,twist_h2,twist_L_hex,twist_is_L_prime,twist_bits_L,
#   method,time_ms

import sys, csv, time, os

inp  = sys.argv[1]
outp = sys.argv[2]

PRP         = (os.environ.get("PRP", "1") == "1")
SKIP_TWIST  = (os.environ.get("SKIP_TWIST", "0") == "1")
SMALL_BOUND = Integer(os.environ.get("SMALL_BOUND", "100000"))

def parse_hex_to_int(hs):
    s = str(hs).strip().lower()
    if s.startswith("0x"):
        return Integer(int(s, 16))
    return Integer(int(s, 16))

def edwards_to_weierstrass(F, a, d):
    r"""
    build Weierstrass model E(Fp) which is isomorphic (over Fp) to the twisted Edwards curve
    """

    # twisted Edwards -> Montgomery
    denom = a - d
    if denom == 0:
        raise ValueError("a == d (singular Edwards)")

    A = F(2) * (a + d) / denom
    B = F(4) / denom
    if B == 0:
        raise ValueError("B == 0 in Edwards->Montgomery")

    # Montgomery: B u^2 = u^3 + A u^2 + u
    # Montgomery -> short Weierstrass
    s = B
    alpha = A / (F(3) * B)
    a4 = s**(-2) - F(3) * alpha**2
    a6 = -alpha * (alpha**2 + a4)

    # returns y^2 = x^3 + a4 x + a6
    E = EllipticCurve(F, [0, 0, 0, a4, a6])
    return E, "edwards-weierstrass"

def build_curve(F, a, d):
    return edwards_to_weierstrass(F, a, d)

def largest_2_part(n):
    """
    split n = 2^v2 * L, returning (v2, h2, L) where h2 = 2^v2 is 2-adic cofactor.
    """
    v2 = valuation(n, 2)
    h2 = Integer(1) << v2
    L  = n >> v2
    return (v2, h2, L)

def strip_small_odd_factors(n, bound):
    """
    remove odd prime factors of L up to this bound
    """
    m = Integer(n)
    if bound <= 2:
        return m
    for q in prime_range(3, int(bound) + 1):
        while m % q == 0:
            m //= q
    return m

def is_prime_fast(n):
    """
    fast primality check for L:
      - is_pseudoprime() when PRP=1
      - is_prime(proof=True) when PRP=0
    """
    if PRP:
        return Integer(n).is_pseudoprime()
    else:
        return Integer(n).is_prime(proof=True)

with open(inp, "r", newline="") as fin, open(outp, "w", newline="") as fout:
    r = csv.DictReader(fin)

    fieldnames = [
        "id","p_hex","a_hex","d_hex",
        "order_hex","h2","L_hex","is_L_prime","bits_L",
        "twist_order_hex","twist_h2","twist_L_hex","twist_is_L_prime","twist_bits_L",
        "method","time_ms"
    ]
    w = csv.DictWriter(fout, fieldnames=fieldnames)
    w.writeheader()

    total = sum(1 for _ in open(inp, "r"))
    fin.seek(0)
    r = csv.DictReader(fin)

    idx = 0
    for row in r:
        idx += 1
        rid   = row.get("id")
        if not rid:
            # if no id column, just use (row_index - 1)
            rid = str(idx - 1)

        p_hex = (row.get("p_hex") or "").strip()
        a_hex = (row.get("a_hex") or "").strip()
        d_hex = (row.get("d_hex") or "").strip()

        t0 = time.time()
        try:
            p = parse_hex_to_int(p_hex)
            F = GF(p)

            a = F(parse_hex_to_int(a_hex))
            d = F(parse_hex_to_int(d_hex))

            # sanity checks - non-singular twisted Edwards
            if a == 0 or d == 0 or a == d:
                raise ValueError("singular candidate (a==0 or d==0 or a==d)")

            # build curve and count #E(Fp)
            E, method = build_curve(F, a, d)
            orderE = E.cardinality()

            # split 2-adic part
            v2, h2, L = largest_2_part(orderE)

            # strip small odd primes from L 
            L = strip_small_odd_factors(L, SMALL_BOUND)

            is_L_prime = 1 if is_prime_fast(L) else 0
            bits_L     = Integer(L).nbits() if L > 0 else 0

            # quadratic twist cardinality
            if SKIP_TWIST:
                twist_order_hex = "0x0"
                twist_h2        = "0"
                twist_L_hex     = "0x0"
                twist_is_prime  = 0
                twist_bits_L    = 0
            else:
                Etw    = E.quadratic_twist()
                orderT = Etw.cardinality()
                tv2, th2, TL = largest_2_part(orderT)
                TL = strip_small_odd_factors(TL, SMALL_BOUND)
                t_is_prime    = 1 if is_prime_fast(TL) else 0

                twist_order_hex = hex(int(orderT))
                twist_h2        = str(int(th2))
                twist_L_hex     = hex(int(TL))
                twist_is_prime  = t_is_prime
                twist_bits_L    = Integer(TL).nbits() if TL > 0 else 0

            w.writerow(dict(
                id = rid,
                p_hex = p_hex,
                a_hex = a_hex,
                d_hex = d_hex,

                order_hex = hex(int(orderE)),
                h2 = str(int(h2)),
                L_hex = hex(int(L)),
                is_L_prime = int(is_L_prime),
                bits_L = int(bits_L),

                twist_order_hex = twist_order_hex,
                twist_h2 = twist_h2,
                twist_L_hex = twist_L_hex,
                twist_is_L_prime = int(twist_is_prime),
                twist_bits_L = int(twist_bits_L),

                method = method,
                time_ms = int(round(1000*(time.time() - t0))),
            ))

        except Exception as e:
            w.writerow(dict(
                id = rid,
                p_hex = p_hex,
                a_hex = a_hex,
                d_hex = d_hex,
                order_hex = "0x0",
                h2 = "0",
                L_hex = "0x0",
                is_L_prime = 0,
                bits_L = 0,
                twist_order_hex = "0x0",
                twist_h2 = "0",
                twist_L_hex = "0x0",
                twist_is_L_prime = 0,
                twist_bits_L = 0,
                method = "error",
                time_ms = int(round(1000*(time.time() - t0))),
            ))