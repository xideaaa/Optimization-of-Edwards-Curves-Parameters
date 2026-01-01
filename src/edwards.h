#ifndef EDWARDS_H
#define EDWARDS_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include "finite_field.h"

// twisted Edwards over Fp: a*x^2 + y^2 = 1 + d*x^2*y^2
// Extended coordinates: P = (X:Y:Z:T),  x = X/Z,  y = Y/Z,  T = (X*Y)/Z

typedef struct {
    const ff_ctx_t* f; // base field 
    ff_elem a; // parameter a
    ff_elem d; // parameter d 
} ed_curve_t;

typedef struct {
    ff_elem X, Y, Z, T; // extended coordinates
} ed_point_t;

bool ed_init(ed_curve_t* ec, const ff_ctx_t* f, const ff_elem* a, const ff_elem* d);
bool ed_init_hex(ed_curve_t* ec, const ff_ctx_t* f, const char* a_hex, const char* d_hex);

// id: (0,1) -> (0:1:1:0)
ed_point_t ed_identity(const ed_curve_t* ec);

// affine validation
bool ed_oncurve_affine(const ff_elem* x, const ff_elem* y, const ed_curve_t* ec);
bool ed_from_affine(ed_point_t* P, const ff_elem* x, const ff_elem* y, const ed_curve_t* ec);
bool ed_to_affine(ff_elem* x, ff_elem* y, const ed_point_t* P, const ed_curve_t* ec);

// group ADD, DOUBLE
void ed_add   (ed_point_t* R, const ed_point_t* P, const ed_point_t* Q, const ed_curve_t* ec);
void ed_double(ed_point_t* R, const ed_point_t* P, const ed_curve_t* ec);

// group SMUL
void ed_scalar_mul_bits(ed_point_t* R,
                        const ed_point_t* P,
                        const uint8_t* k_le_bits, size_t nbits,
                        const ed_curve_t* ec);
void ed_scalar_mul_u64(ed_point_t* R, const ed_point_t* P, uint64_t k, const ed_curve_t* ec);

void ed_scalar_mul(ed_point_t* R, const ed_point_t* P, uint64_t k, const ed_curve_t* ec);

typedef enum {
    ED_SMUL_BINARY = 0,  
    ED_SMUL_WNAF   = 1,   
    ED_SMUL_FIXED  = 2    
} ed_smul_algo_t;

typedef struct {
    int w; // window-size
    size_t ntable;       
    ed_point_t* tbl; // basis {1*G, 2*G, ..., (2^w-1)*G}
} ed_fixed_precomp_t;


void ed_scalar_mul_wnaf_bits(ed_point_t* R,
                             const ed_point_t* P,
                             const uint8_t* k_le_bits, size_t nbits,
                             int w,
                             const ed_curve_t* ec);

int  ed_fixed_precompute(ed_fixed_precomp_t* pc,
                         const ed_point_t* G, int w,
                         const ed_curve_t* ec);

void ed_scalar_mul_fixed_bits(ed_point_t* R,
                              const uint8_t* k_le_bits, size_t nbits,
                              const ed_fixed_precomp_t* pc,
                              const ed_curve_t* ec);

void ed_fixed_precomp_free(ed_fixed_precomp_t* pc);

void ed_scalar_mul_ex(ed_point_t* R,
                      const ed_point_t* P,
                      const uint8_t* k_le_bits, size_t nbits,
                      ed_smul_algo_t algo,
                      const ed_fixed_precomp_t* pc,
                      const ed_curve_t* ec);

#endif /* EDWARDS_H */
