#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// fixed multi-precision 8 limbs * 64 = 512 bits
#define FF_NLIMBS 8

// Fp element (little-endian limbs: v[0] LSB)
typedef struct { uint64_t v[FF_NLIMBS]; } ff_elem;

// field context: modulus p, montgomery constants, limb count
typedef struct {
    uint64_t p[FF_NLIMBS]; // modulus
    uint64_t n0; // -p[0]^{-1} mod 2^64
    uint64_t R2[FF_NLIMBS]; // R^2 mod p; R = 2^(64*FF_NLIMBS)
    int      nlimbs; // limb count of p (1,...,8)
} ff_ctx_t;


bool ff_init(ff_ctx_t* ctx, uint64_t p_u64);

bool ff_init_hex(ff_ctx_t* ctx, const char* hex_p);

static inline const uint64_t* ff_modulus_limbs(const ff_ctx_t* ctx) { return ctx->p; }
static inline int ff_modulus_nlimbs(const ff_ctx_t* ctx) { return ctx->nlimbs; }

static inline uint64_t ff_modulus(const ff_ctx_t* ctx) { return ctx->p[0]; }

bool ff_from_hex(const char* hex, ff_elem* a);
void ff_to_hex(const ff_elem* a, char* out, size_t out_sz);
void ff_from_u64(uint64_t x, ff_elem* a);
bool ff_to_u64(const ff_elem* a, uint64_t* out);

// 512-bit element arithmetic
void ff_red_e(ff_elem* a, const ff_ctx_t* ctx);
void ff_add_e(ff_elem* r, const ff_elem* a, const ff_elem* b, const ff_ctx_t*);
void ff_sub_e(ff_elem* r, const ff_elem* a, const ff_elem* b, const ff_ctx_t*);
void ff_neg_e(ff_elem* r, const ff_elem* a, const ff_ctx_t*);
void ff_mul_e(ff_elem* r, const ff_elem* a, const ff_elem* b, const ff_ctx_t*);
void ff_pow_e(ff_elem* r, const ff_elem* a, const ff_elem* e, const ff_ctx_t*);
bool ff_inv_e(ff_elem* r, const ff_elem* a, const ff_ctx_t*);

// legacy 64-bit API
uint64_t ff_red(uint64_t a, const ff_ctx_t* ctx);
uint64_t ff_add(uint64_t a, uint64_t b, const ff_ctx_t* ctx);
uint64_t ff_sub(uint64_t a, uint64_t b, const ff_ctx_t* ctx);
uint64_t ff_neg(uint64_t a, const ff_ctx_t* ctx);
uint64_t ff_mul(uint64_t a, uint64_t b, const ff_ctx_t* ctx);
uint64_t ff_pow(uint64_t a, uint64_t e, const ff_ctx_t* ctx);
uint64_t ff_inv(uint64_t a, const ff_ctx_t* ctx);

static inline int ff_eq(uint64_t a, uint64_t b, const ff_ctx_t* ctx) {
    (void)ctx;
    uint64_t x = a ^ b;
    x |= (uint64_t)-(int64_t)x;
    x >>= 63;
    return (int)(1 ^ x);
}

#ifdef __cplusplus
}
#endif

#endif /* FINITE_FIELD_H */
