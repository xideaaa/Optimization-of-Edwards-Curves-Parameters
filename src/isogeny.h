#ifndef ISOGENY_H
#define ISOGENY_H

#include "finite_field.h"
#include "cost.h"

typedef struct {
    ff_elem a;
    ff_elem d;
} ed_curve_params_t;

void ed_curve_params_from_hex(ed_curve_params_t* C,
                              const char* a_hex,
                              const char* d_hex);

void ed_iso4_curve(ed_curve_params_t* out,
                   const ed_curve_params_t* in,
                   const ff_ctx_t* F);

void ed_iso4_map(ff_elem* x2, ff_elem* y2,
                 const ff_elem* x,
                 const ff_elem* y,
                 const ed_curve_params_t* C,
                 const ff_ctx_t* F);

#endif
