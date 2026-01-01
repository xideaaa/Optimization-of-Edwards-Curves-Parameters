#ifndef ED_UTILS_H
#define ED_UTILS_H

#include <stdint.h>
#include "finite_field.h"
#include "edwards.h"

static inline int ff_eq_zero(const ff_elem* x){
    for (int i=0;i<FF_NLIMBS;i++) if (x->v[i]) return 0;
    return 1;
}
static inline int ff_eq_one(const ff_elem* x){
    if (x->v[0] != 1) return 0;
    for (int i=1;i<FF_NLIMBS;i++) if (x->v[i]) return 0;
    return 1;
}

typedef struct { uint64_t s; } ed_rng_t;
static inline void ed_rng_seed(ed_rng_t* r, uint64_t x){
    r->s = x ? x : 0x9e3779b97f4a7c15ull;
}
static inline uint64_t ed_rng_next(ed_rng_t* r){
    uint64_t x = r->s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    return r->s = x;
}

// Legendre symbol (a|p): 1 residue, 0 zero, -1 non-residue
int ff_legendre(const ff_elem* a, const ff_ctx_t* F);

static inline int ff_is_square(const ff_elem* a, const ff_ctx_t* F){
    int L = ff_legendre(a, F);
    return (L == 1);
}

// a!=0, d!=0, a!=d, a is square, d is non-square
int ed_complete_candidate(const ff_elem* a, const ff_elem* d, const ff_ctx_t* F);

void ff_random(ff_elem* r, ed_rng_t* R, const ff_ctx_t* F);

// Tonelli–Shanks square root: returns 1 and sets r if sqrt exists, 0 otherwise
int ff_sqrt_ts(ff_elem* r, const ff_elem* a, const ff_ctx_t* F);

#endif /* ED_UTILS_H */
