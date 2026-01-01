#include "edwards.h"
#include "cost.h"
#include <string.h>

static void ff_set0(ff_elem* r){ memset(r->v, 0, sizeof r->v); }
static void ff_set1(ff_elem* r){ memset(r->v, 0, sizeof r->v); r->v[0] = 1; }
static void ff_set2(ff_elem* r){ memset(r->v, 0, sizeof r->v); r->v[0] = 2; }

bool ed_init(ed_curve_t* ec, const ff_ctx_t* f, const ff_elem* a, const ff_elem* d){
    if (!ec || !f || !a || !d) return false;
    ec->f = f;
    ec->a = *a;
    ec->d = *d;

    ff_elem z; ff_set0(&z);
    if (memcmp(ec->a.v, z.v, sizeof z.v) == 0) return false;
    if (memcmp(ec->d.v, z.v, sizeof z.v) == 0) return false;
    if (memcmp(ec->a.v, ec->d.v, sizeof z.v) == 0) return false;
    return true;
}

bool ed_init_hex(ed_curve_t* ec, const ff_ctx_t* f, const char* a_hex, const char* d_hex){
    if (!ec || !f || !a_hex || !d_hex) return false;
    ff_elem a,d;
    if (!ff_from_hex(a_hex,&a)) return false;
    if (!ff_from_hex(d_hex,&d)) return false;
    return ed_init(ec,f,&a,&d);
}

ed_point_t ed_identity(const ed_curve_t* ec){
    (void)ec;
    ed_point_t I;
    ff_set0(&I.X);
    ff_set1(&I.Y);
    ff_set1(&I.Z);
    ff_set0(&I.T);
    return I;
}

bool ed_oncurve_affine(const ff_elem* x, const ff_elem* y, const ed_curve_t* ec){
    const ff_ctx_t* Fctx = ec->f;
    ff_elem x2,y2,lhs,rhs,one,tmp;

    COUNTED_MUL(&x2,x,x,Fctx); // x^2
    COUNTED_MUL(&y2,y,y,Fctx); //* y^2 

    // lhs = a*x^2 + y^2
    COUNTED_MUL(&lhs,&ec->a,&x2,Fctx);
    COUNTED_ADD(&lhs,&lhs,&y2,Fctx);

    // rhs = 1 + d*x^2*y^2
    ff_set1(&one);
    COUNTED_MUL(&tmp,&x2,&y2,Fctx); // x^2*y^2 
    COUNTED_MUL(&rhs,&ec->d,&tmp,Fctx); // d*x^2*y^2 
    COUNTED_ADD(&rhs,&one,&rhs,Fctx); // 1 + d*x^2*y^2

    return memcmp(lhs.v, rhs.v, sizeof lhs.v) == 0;
}

// affine <-> extended
bool ed_from_affine(ed_point_t* P, const ff_elem* x, const ff_elem* y, const ed_curve_t* ec){
    const ff_ctx_t* Fctx = ec->f;
    (void)Fctx;
    P->X = *x;
    P->Y = *y;
    ff_set1(&P->Z);
    COUNTED_MUL(&P->T, x, y, Fctx); // T = x*y
    return true;
}

bool ed_to_affine(ff_elem* x, ff_elem* y, const ed_point_t* P, const ed_curve_t* ec){
    const ff_ctx_t* Fctx = ec->f;
    ff_elem zinv;
    if (!ff_inv_e(&zinv, &P->Z, Fctx))
        return false; // Z == 0 (point at infinity)

    COUNTED_MUL(x, &P->X, &zinv, Fctx);
    COUNTED_MUL(y, &P->Y, &zinv, Fctx);
    return true;
}

/* Unified Addition
 *
 *  A = X1*X2
 *  B = Y1*Y2
 *  C = d*T1*T2
 *  D = Z1*Z2
 *  E = (X1+Y1)(X2+Y2) - A - B
 *  F = D - C
 *  G = D + C
 *  H = B - a*A
 *  X3 = E*F
 *  Y3 = G*H
 *  T3 = E*H
 *  Z3 = F*G
 */
void ed_add(ed_point_t* R,
            const ed_point_t* P,
            const ed_point_t* Q,
            const ed_curve_t* ec)
{
    const ff_ctx_t* Fctx = ec->f;
    ff_elem A,B,C,D,E,Fm,G,H;
    ff_elem X1pY1, X2pY2;

    // A = X1*X2
    COUNTED_MUL(&A, &P->X, &Q->X, Fctx);

    // B = Y1*Y2
    COUNTED_MUL(&B, &P->Y, &Q->Y, Fctx);

    // C = d*T1*T2 
    COUNTED_MUL(&C, &P->T, &Q->T, Fctx);
    COUNTED_MUL(&C, &C, &ec->d, Fctx);

    // D = Z1*Z2
    COUNTED_MUL(&D, &P->Z, &Q->Z, Fctx);

    // E = (X1+Y1)(X2+Y2) - A - B 
    COUNTED_ADD(&X1pY1, &P->X, &P->Y, Fctx);
    COUNTED_ADD(&X2pY2, &Q->X, &Q->Y, Fctx);
    COUNTED_MUL(&E, &X1pY1, &X2pY2, Fctx);
    COUNTED_SUB(&E, &E, &A, Fctx);
    COUNTED_SUB(&E, &E, &B, Fctx);

    // F = D - C
    COUNTED_SUB(&Fm, &D, &C, Fctx);

    // G = D + C 
    COUNTED_ADD(&G, &D, &C, Fctx);

    // H = B - a*A 
    COUNTED_MUL(&H, &ec->a, &A, Fctx);
    COUNTED_SUB(&H, &B, &H, Fctx);

    // X3 = E*F 
    COUNTED_MUL(&R->X, &E, &Fm, Fctx);

    // Y3 = G*H 
    COUNTED_MUL(&R->Y, &G, &H, Fctx);

    // T3 = E*H 
    COUNTED_MUL(&R->T, &E, &H, Fctx);

    // Z3 = F*G 
    COUNTED_MUL(&R->Z, &Fm, &G, Fctx);
}

/* Doubling 
 *
 *  A = X^2
 *  B = Y^2
 *  C = d*T^2
 *  D = Z^2
 *  E = (X+Y)^2 - A - B
 *  F = D - C
 *  G = D + C
 *  H = B - a*A
 *  X3 = E*F
 *  Y3 = G*H
 *  T3 = E*H
 *  Z3 = F*G
 */
void ed_double(ed_point_t* R,
               const ed_point_t* P,
               const ed_curve_t* ec)
{
    const ff_ctx_t* Fctx = ec->f;
    ff_elem A,B,C,D,E,Fm,G,H,XpY;

    // A = X^2 
    COUNTED_MUL(&A, &P->X, &P->X, Fctx);

    // B = Y^2 
    COUNTED_MUL(&B, &P->Y, &P->Y, Fctx);

    // C = d*T^2 
    COUNTED_MUL(&C, &P->T, &P->T, Fctx);
    COUNTED_MUL(&C, &C, &ec->d, Fctx);

    // D = Z^2 
    COUNTED_MUL(&D, &P->Z, &P->Z, Fctx);

    // E = (X+Y)^2 - A - B 
    COUNTED_ADD(&XpY, &P->X, &P->Y, Fctx);
    COUNTED_MUL(&E, &XpY, &XpY, Fctx);
    COUNTED_SUB(&E, &E, &A, Fctx);
    COUNTED_SUB(&E, &E, &B, Fctx);

    // F = D - C 
    COUNTED_SUB(&Fm, &D, &C, Fctx);

    // G = D + C 
    COUNTED_ADD(&G, &D, &C, Fctx);

    // H = B - a*A 
    COUNTED_MUL(&H, &ec->a, &A, Fctx);
    COUNTED_SUB(&H, &B, &H, Fctx);

    // X3 = E*F 
    COUNTED_MUL(&R->X, &E, &Fm, Fctx);

    // Y3 = G*H 
    COUNTED_MUL(&R->Y, &G, &H, Fctx);

    // T3 = E*H 
    COUNTED_MUL(&R->T, &E, &H, Fctx);

    // Z3 = F*G
    COUNTED_MUL(&R->Z, &Fm, &G, Fctx);
}


static void cswap(ed_point_t* A, ed_point_t* B, uint64_t mask){
    for (int i = 0; i < FF_NLIMBS; i++){
#define SWAP_FIELD(F) do{ \
    uint64_t t = (A->F.v[i] ^ B->F.v[i]) & mask; \
    A->F.v[i] ^= t; \
    B->F.v[i] ^= t; \
}while(0)
        SWAP_FIELD(X);
        SWAP_FIELD(Y);
        SWAP_FIELD(Z);
        SWAP_FIELD(T);
#undef SWAP_FIELD
    }
}

// SMUL - Montgomery ladder
void ed_scalar_mul_bits(ed_point_t* R,
                        const ed_point_t* P,
                        const uint8_t* k_le_bits, size_t nbits,
                        const ed_curve_t* ec)
{
    ed_point_t R0 = ed_identity(ec);
    ed_point_t R1 = *P;

    if (nbits == 0) {
        *R = R0;
        return;
    }

    for (ssize_t i = (ssize_t)nbits - 1; i >= 0; --i) {
        uint8_t byte = k_le_bits[i >> 3];
        uint64_t b = (byte >> (i & 7)) & 1u;
        uint64_t mask = (uint64_t)-(b & 1u);

        cswap(&R0, &R1, mask);
        ed_add(&R1, &R0, &R1, ec);
        ed_double(&R0, &R0, ec);
        cswap(&R0, &R1, mask);
    }
    *R = R0;
}

void ed_scalar_mul_u64(ed_point_t* R, const ed_point_t* P, uint64_t k, const ed_curve_t* ec){
    uint8_t bits[8] = {0};
    for (int i = 0; i < 64; i++){
        if ((k >> i) & 1u)
            bits[i >> 3] |= (uint8_t)(1u << (i & 7));
    }
    ed_scalar_mul_bits(R, P, bits, 64, ec);
}

void ed_scalar_mul(ed_point_t* R, const ed_point_t* P, uint64_t k, const ed_curve_t* ec){
    ed_scalar_mul_u64(R, P, k, ec);
}
