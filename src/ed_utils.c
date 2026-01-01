#include <string.h>
#include "ed_utils.h"

// compute Legendre symbol using Eulers criterion
int ff_legendre(const ff_elem* a, const ff_ctx_t* F){
    if (ff_eq_zero(a)) return 0;

    // e = (p-1)/2
    ff_elem e;
    memcpy(e.v, F->p, sizeof e.v);

    // e = p-1
    uint64_t borrow = 0;
    e.v[0] = (uint64_t)((__uint128_t)e.v[0] - 1 - borrow);
    borrow = (e.v[0] > F->p[0]) ? 1 : 0;
    for (int i = 1; i < FF_NLIMBS; ++i) {
        uint64_t before = e.v[i];
        e.v[i] = (uint64_t)((__uint128_t)e.v[i] - 0 - borrow);
        borrow = (before < borrow) ? 1 : 0;
    }
    // e >>= 1 
    for (int i = 0; i < FF_NLIMBS; ++i){
        uint64_t next = (i+1 < FF_NLIMBS) ? (e.v[i+1] & 1u) : 0u;
        e.v[i] = (e.v[i] >> 1) | (next << 63);
    }

    ff_elem r;
    ff_pow_e(&r, a, &e, F); // r = a^((p-1)/2) mod p 

    if (ff_eq_one(&r)) return 1;

    // pm1 = p-1
    ff_elem pm1;
    memcpy(pm1.v, F->p, sizeof pm1.v);
    uint64_t br = 0;
    pm1.v[0] = (uint64_t)((__uint128_t)pm1.v[0] - 1 - br);
    br = (pm1.v[0] > F->p[0]) ? 1 : 0;
    for (int i = 1; i < FF_NLIMBS; ++i) {
        uint64_t before = pm1.v[i];
        pm1.v[i] = (uint64_t)((__uint128_t)pm1.v[i] - 0 - br);
        br = (before < br) ? 1 : 0;
    }

    int is_pm1 = 1;
    for (int i=0;i<FF_NLIMBS;i++) if (r.v[i] != pm1.v[i]) { is_pm1 = 0; break; }
    return is_pm1 ? -1 : -1;
}

int ed_complete_candidate(const ff_elem* a, const ff_elem* d, const ff_ctx_t* F){
    if (ff_eq_zero(a) || ff_eq_zero(d)) return 0;
    int eq = 1; for (int i=0;i<FF_NLIMBS;i++) if (a->v[i] != d->v[i]) { eq = 0; break; }
    if (eq) return 0;
    if (!ff_is_square(a,F)) return 0;
    if (ff_is_square(d,F))  return 0;
    return 1;
}

static inline int ff_is_odd(const ff_elem* a){ return (a->v[0] & 1u) != 0; }

static void ff_rshift1(ff_elem* x){
    for (int i = 0; i < FF_NLIMBS; ++i){
        uint64_t next = (i+1 < FF_NLIMBS) ? (x->v[i+1] & 1u) : 0u;
        x->v[i] = (x->v[i] >> 1) | (next << 63);
    }
}

static void ff_copy(ff_elem* dst, const ff_elem* src){ memcpy(dst->v, src->v, sizeof dst->v); }

void ff_random(ff_elem* r, ed_rng_t* R, const ff_ctx_t* F){
    for (int i=0;i<FF_NLIMBS;i++) r->v[i] = ed_rng_next(R);
    ff_red_e(r, F);
}

// Tonelli–Shanks (TS) square root in Fp
int ff_sqrt_ts(ff_elem* out, const ff_elem* a, const ff_ctx_t* F){
    if (ff_eq_zero(a)){ memset(out->v,0,sizeof out->v); return 1; }
    if (!ff_is_square(a, F)) return 0;

    // p-1 = q * 2^s, q odd
    ff_elem pm1; memcpy(pm1.v, F->p, sizeof pm1.v);
    uint64_t br=0;
    pm1.v[0] = (uint64_t)((__uint128_t)pm1.v[0] - 1 - br);
    br = (pm1.v[0] > F->p[0]) ? 1 : 0;
    for (int i=1;i<FF_NLIMBS;i++){ uint64_t t=pm1.v[i]; pm1.v[i] = (uint64_t)((__uint128_t)pm1.v[i] - 0 - br); br=(t<br); }

    ff_elem q; ff_copy(&q,&pm1);
    int s = 0;
    while (!ff_is_odd(&q)){ ff_rshift1(&q); s++; }

    // find non-residue z 
    ed_rng_t RR; ed_rng_seed(&RR, 0x123456789abcdefull ^ F->p[0]);
    ff_elem z;
    for (;;){
        ff_random(&z, &RR, F);
        int leg = ff_legendre(&z,F);
        if (leg == -1) break;
    }

    // TS alg
    ff_elem c, t, r, b, tmp;
    ff_pow_e(&c, &z, &q, F); //  c = z^q 

    // r = a^((q+1)/2)
    ff_elem qp1; ff_copy(&qp1, &q);
    uint64_t carry = 1; // +1
    for (int i=0;i<FF_NLIMBS;i++){
        __uint128_t ssum = (__uint128_t)qp1.v[i] + carry;
        qp1.v[i] = (uint64_t)ssum;
        carry = (uint64_t)(ssum >> 64);
    }
    ff_rshift1(&qp1);
    ff_pow_e(&r, a, &qp1, F);

    ff_pow_e(&t, a, &q, F);
    int m = s;

    ff_elem one; memset(one.v,0,sizeof one.v); one.v[0]=1;

    while (1){
        if (ff_eq_one(&t)){ ff_copy(out,&r); return 1; }
        // find smallest i in [1,m) s.t. t^(2^i) == 1
        int i=1;
        ff_copy(&tmp, &t);
        for (; i<m; ++i){
            ff_mul_e(&tmp, &tmp, &tmp, F);
            if (ff_eq_one(&tmp)) break;
        }
        // b = c^(2^(m-i-1))
        int e = m - i - 1;
        ff_copy(&b, &c);
        for (int j=0;j<e;j++){
            ff_mul_e(&b, &b, &b, F);
        }
        // r = r*b; c = b^2; t = t*c
        ff_mul_e(&r, &r, &b, F);
        ff_mul_e(&b, &b, &b, F);
        ff_copy(&c, &b);
        ff_mul_e(&t, &t, &c, F);
        m = i;
    }
}
