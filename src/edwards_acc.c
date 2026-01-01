#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "edwards_acc.h"  
#include "cost.h"        

static inline uint64_t ct_eq_u64_mask(uint64_t a, uint64_t b){
    uint64_t x = a ^ b;
    x |= (0u - x);
    return ~(x >> 63);
}

static inline uint64_t ct_mask_from_bool(int cond){
    return (uint64_t)(-(cond != 0));
}

static inline void ed_ct_select(ed_point_t* R, const ed_point_t* A, uint64_t mask){
    for (int i=0;i<FF_NLIMBS;i++){
        R->X.v[i] = (R->X.v[i] & ~mask) | (A->X.v[i] & mask);
        R->Y.v[i] = (R->Y.v[i] & ~mask) | (A->Y.v[i] & mask);
        R->Z.v[i] = (R->Z.v[i] & ~mask) | (A->Z.v[i] & mask);
        R->T.v[i] = (R->T.v[i] & ~mask) | (A->T.v[i] & mask);
    }
}

// -P = (-X, Y, Z, -T)
static inline void ed_neg_local(ed_point_t* R, const ed_point_t* P, const ed_curve_t* ec){
    const ff_ctx_t* Fctx = ec->f;              
    COUNTED_NEG(&R->X, &P->X, Fctx);
    COUNTED_NEG(&R->T, &P->T, Fctx);
    R->Y = P->Y;
    R->Z = P->Z;
}

static inline void ed_set_identity(ed_point_t* R, const ed_curve_t* ec){
    ff_elem x0, y1;
    ff_from_u64(0, &x0);
    ff_from_u64(1, &y1);
    ed_from_affine(R, &x0, &y1, ec);
}

static inline int get_bit_lsb_first(const uint8_t* bits, size_t i){
    return (bits[i>>3] >> (i & 7)) & 1;
}

// w-NAF
typedef struct {
    int8_t  digits[520]; // basis {0, =/-1, +/-3, ..., +/-(2^{w-1}-1)}
    size_t  len;
    int     w;
} ed_wnaf_t;

static void ed_wnaf_encode(ed_wnaf_t* out, const uint8_t* k_bits, size_t nbits, int w){
    assert(out && k_bits);
    assert(w >= 2 && w <= 7);

    const size_t nbytes = (nbits + 7) >> 3;
    uint8_t* k = (uint8_t*)malloc(nbytes ? nbytes : 1);
    if (!k) { out->len = 0; out->w = w; return; }
    memcpy(k, k_bits, nbytes);

    size_t i = 0, out_len = 0;
    const uint32_t two_w   = 1u << w;
    const uint32_t two_wm1 = 1u << (w - 1);

    while (i < nbits) {
        int di = 0;
        if (get_bit_lsb_first(k, i) == 1){
            uint32_t limb = 0;
            size_t upto = (i + (size_t)w <= nbits) ? (size_t)w : (nbits - i);
            for (size_t b = 0; b < upto; ++b){
                limb |= (uint32_t)get_bit_lsb_first(k, i + b) << b;
            }
            int u = (int)(limb & (two_w - 1));
            if (u >= (int)two_wm1) u -= (int)two_w;
            if ((u & 1) == 0) { u += (u > 0) ? -1 : 1; }
            di = u;

            // k = k - u * 2^i
            int absu = (u < 0) ? -u : u;
            for (int b = 0; b < w; ++b){
                if (((absu >> b) & 1) == 0) continue;
                size_t pos = i + (size_t)b;
                int carry = (u > 0) ? -1 : +1;
                while (carry != 0 && pos < nbits){
                    int bit = get_bit_lsb_first(k, pos);
                    int sum = bit + carry;
                    int newbit = sum & 1;
                    carry = (sum >> 1);
                    if (newbit) k[pos>>3] |=  (1u << (pos & 7));
                    else        k[pos>>3] &= ~(1u << (pos & 7));
                    pos++;
                }
            }
        } else {
            di = 0;
        }

        out->digits[i] = (int8_t)di;
        out_len = i + 1;

        // k >>= 1 (bit)
        uint8_t carry = 0;
        for (size_t b = (nbytes ? nbytes : 1); b-- > 0; ){
            uint8_t nb = (uint8_t)((k[b] >> 1) | (carry << 7));
            carry = (uint8_t)(k[b] & 1);
            k[b] = nb;
        }
        i++;
    }

    free(k);
    out->len = out_len;
    out->w   = w;
}

void ed_scalar_mul_wnaf_bits(ed_point_t* R,
                             const ed_point_t* P,
                             const uint8_t* k_le_bits, size_t nbits,
                             int w,
                             const ed_curve_t* ec)
{
    assert(R && P && k_le_bits && ec);
    assert(w >= 2 && w <= 7);

    ed_wnaf_t naf;
    ed_wnaf_encode(&naf, k_le_bits, nbits, w);

    // build basis T[i] = (2*i+1)P, i = 0..2^{w-1}-1
    const int max_odd = (1 << (w-1)) - 1;
    const int tbl_size = (max_odd > 0) ? max_odd : 1;
    ed_point_t* T = (ed_point_t*)malloc(sizeof(ed_point_t) * (size_t)tbl_size);
    if (!T) { ed_scalar_mul_bits(R, P, k_le_bits, nbits, ec); return; }

    ed_point_t P2;
    ed_double(&P2, P, ec);
    T[0] = *P;
    for (int i = 1; i < tbl_size; ++i){
        ed_add(&T[i], &T[i-1], &P2, ec);
    }

    ed_set_identity(R, ec);

    for (size_t i = naf.len; i-- > 0; ){
        ed_double(R, R, ec);

        int8_t di = naf.digits[i];
        int sign = (di < 0);
        int absd = sign ? -di : di;

        if (absd != 0){
            int idx = (absd - 1) >> 1;

            ed_point_t pick = T[0];
            for (int j = 0; j < tbl_size; ++j){
                uint64_t m = ct_eq_u64_mask((uint64_t)j, (uint64_t)idx);
                ed_ct_select(&pick, &T[j], m);
            }

            ed_point_t negpick;
            ed_neg_local(&negpick, &pick, ec);
            ed_point_t addend = pick;
            uint64_t mneg = ct_mask_from_bool(sign);
            ed_ct_select(&addend, &negpick, mneg);

            ed_add(R, R, &addend, ec);
        }
    }

    free(T);
}

// fixed-window
int ed_fixed_precompute(ed_fixed_precomp_t* pc,
                        const ed_point_t* G, int w,
                        const ed_curve_t* ec)
{
    if (!pc || !G || !ec) return 0;
    if (w < 2 || w > 7)   return 0;

    pc->w = w;
    pc->ntable = ((size_t)1u << w) - 1;       
    pc->tbl = (ed_point_t*)malloc(sizeof(ed_point_t) * pc->ntable);
    if (!pc->tbl) { pc->ntable = 0; return 0; }

    // basis  {1*G, 2*G, 3*G, ..., (2^w-1)*G}
    pc->tbl[0] = *G;
    for (size_t i = 1; i < pc->ntable; ++i){
        ed_add(&pc->tbl[i], &pc->tbl[i-1], G, ec);
    }
    return 1;
}

void ed_fixed_precomp_free(ed_fixed_precomp_t* pc){
    if (!pc) return;
    if (pc->tbl){
        memset(pc->tbl, 0, sizeof(ed_point_t)*pc->ntable);
        free(pc->tbl);
    }
    pc->tbl = NULL;
    pc->ntable = 0;
    pc->w = 0;
}

void ed_scalar_mul_fixed_bits(ed_point_t* R,
                              const uint8_t* k_le_bits, size_t nbits,
                              const ed_fixed_precomp_t* pc,
                              const ed_curve_t* ec)
{
    assert(R && pc && pc->tbl && k_le_bits && ec);
    const int w = pc->w;
    const size_t ntable = pc->ntable;

    ed_set_identity(R, ec);

    size_t i = nbits;
    while (i > 0){
        // take window [i-w, i) as unsigned value v (LSB-first)
        uint32_t v = 0;
        size_t take = (i >= (size_t)w) ? (size_t)w : i;
        size_t start = i - take;
        for (size_t b = 0; b < take; ++b){
            int bit = get_bit_lsb_first(k_le_bits, start + b);
            v |= (uint32_t)bit << b;
        }

        // w doubling each step
        for (int d=0; d<w; ++d) ed_double(R, R, ec);

        if (v){
            // lookup of (v * G): v \in [1..2^w-1]
            size_t idx = (size_t)(v - 1);
            ed_point_t pick = pc->tbl[0];
            for (size_t j=0; j<ntable; ++j){
                uint64_t m = ct_eq_u64_mask(j, idx);
                ed_ct_select(&pick, &pc->tbl[j], m);
            }
            ed_add(R, R, &pick, ec);
        }

        i = start;
    }
}

void ed_scalar_mul_ex(ed_point_t* R,
                      const ed_point_t* P,
                      const uint8_t* k_le_bits, size_t nbits,
                      ed_smul_algo_t algo,
                      const ed_fixed_precomp_t* pc,
                      const ed_curve_t* ec)
{
    switch (algo){
        case ED_SMUL_WNAF:
            ed_scalar_mul_wnaf_bits(R, P, k_le_bits, nbits, 5, ec); // default w=5
            break;
        case ED_SMUL_FIXED:
            ed_scalar_mul_fixed_bits(R, k_le_bits, nbits, pc, ec);
            break;
        case ED_SMUL_BINARY:
        default:
            ed_scalar_mul_bits(R, P, k_le_bits, nbits, ec);
            break;
    }
}
