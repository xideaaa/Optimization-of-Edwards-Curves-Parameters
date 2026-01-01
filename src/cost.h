#ifndef COST_H
#define COST_H

#include <stdint.h>
#include <string.h>
#include "finite_field.h"

#ifdef FF_COUNT_OPS

extern _Thread_local uint64_t FF_COUNT_MULS;
extern _Thread_local uint64_t FF_COUNT_SQRS;
extern _Thread_local uint64_t FF_COUNT_ADDS;

static inline void ff_cost_reset(void) {
    FF_COUNT_MULS = 0;
    FF_COUNT_SQRS = 0;
    FF_COUNT_ADDS = 0;
}

typedef struct { uint64_t M, S, A; } ff_cost_t;

static inline ff_cost_t ff_cost_snapshot(void) {
    ff_cost_t c = { FF_COUNT_MULS, FF_COUNT_SQRS, FF_COUNT_ADDS };
    return c;
}

static inline int ff_elem_eq_small(const ff_elem* a, uint64_t small) {
    if (a->v[0] != small) return 0;
    for (int i = 1; i < FF_NLIMBS; ++i) if (a->v[i] != 0) return 0;
    return 1;
}

static inline void COUNTED_MUL(ff_elem* r, const ff_elem* x, const ff_elem* y, const ff_ctx_t* F) {
    if (x == y) {
        ff_mul_e(r, x, y, F);
        FF_COUNT_SQRS++;
    } else {
        ff_mul_e(r, x, y, F);
        FF_COUNT_MULS++;
    }
}

static inline void COUNTED_ADD(ff_elem* r, const ff_elem* x, const ff_elem* y, const ff_ctx_t* F) {
    (void)F;
    ff_add_e(r, x, y, F);
    FF_COUNT_ADDS++;
}

static inline void COUNTED_SUB(ff_elem* r, const ff_elem* x, const ff_elem* y, const ff_ctx_t* F) {
    (void)F;
    ff_sub_e(r, x, y, F);
    FF_COUNT_ADDS++; 
}

static inline void COUNTED_NEG(ff_elem* r, const ff_elem* x, const ff_ctx_t* F) {
    (void)F;
    ff_neg_e(r, x, F);
}

#else

#define ff_cost_reset() ((void)0)
typedef struct { uint64_t M,S,A; } ff_cost_t;
static inline ff_cost_t ff_cost_snapshot(void){ ff_cost_t z={0,0,0}; return z; }

#define COUNTED_MUL  ff_mul_e
#define COUNTED_ADD  ff_add_e
#define COUNTED_SUB  ff_sub_e
#define COUNTED_NEG  ff_neg_e

#endif 

#endif 
