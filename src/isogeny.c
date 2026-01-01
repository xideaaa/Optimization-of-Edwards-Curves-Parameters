#include "isogeny.h"

void ed_curve_params_from_hex(ed_curve_params_t* C,
                              const char* a_hex,
                              const char* d_hex)
{
    ff_from_hex(a_hex, &C->a);
    ff_from_hex(d_hex, &C->d);
}

void ed_iso4_curve(ed_curve_params_t* out,
                   const ed_curve_params_t* in,
                   const ff_ctx_t* F)
{
    // a'= -a
    COUNTED_NEG(&out->a, &in->a, F);

    // d'= d - a
    COUNTED_SUB(&out->d, &in->d, &in->a, F);
}

// f(x,y) = ( 2xy / (y^2 - a x^2), (y^2 + a x^2) / (2 - y^2 - a x^2) )
void ed_iso4_map(ff_elem* x2, ff_elem* y2,
                 const ff_elem* x,
                 const ff_elem* y,
                 const ed_curve_params_t* C,
                 const ff_ctx_t* F)
{
    ff_elem x2_tmp, y2_tmp;
    ff_elem xx, yy, ax2;
    ff_elem num, den;
    ff_elem two;

    ff_from_u64(2, &two);

    // xx = x^2, yy = y^2
    COUNTED_MUL(&xx, x, x, F);
    COUNTED_MUL(&yy, y, y, F);

    // ax2 = a * x^2
    COUNTED_MUL(&ax2, &C->a, &xx, F);

    // x' = 2xy / (y^2 - a x^2)
    COUNTED_MUL(&num, x, y, F);
    COUNTED_MUL(&num, &num, &two, F); // 2xy

    COUNTED_SUB(&den, &yy, &ax2, F); // y^2 - a x^2 
    ff_inv_e(&den, &den, F);

    COUNTED_MUL(&x2_tmp, &num, &den, F);

    // y' = (y^2 + a x^2) / (2 - y^2 - a x^2)
    COUNTED_ADD(&num, &yy, &ax2, F); // y^2 + a x^2

    COUNTED_SUB(&den, &two, &yy, F); // 2 - y^2 
    COUNTED_SUB(&den, &den, &ax2, F); // 2 - y^2 - a x^2
    ff_inv_e(&den, &den, F);

    COUNTED_MUL(&y2_tmp, &num, &den, F);

    *x2 = x2_tmp;
    *y2 = y2_tmp;
}
