// curve_select.c
//
// Selector for twisted Edwards curves
// Input:  counted.csv  (from count_with_sage.sage, SMALL_BOUND=2)
//
//  Header (counted.csv):
//    id,p_hex,a_hex,d_hex,order_hex,h2,L_hex,is_L_prime,bits_L,
//    twist_order_hex,twist_h2,twist_L_hex,twist_is_L_prime,twist_bits_L,
//    method,time_ms
//
// Select curves with:
//   - is_L_prime = 1
//   - bits_L >= BITS_MIN
//   - h2 \in {1,4,8}      
//
// For each selected curve:
//   - init field F_p
//   - init twisted Edwards curve
//   - sample random points P on Edwards
//   - compute Q = h2 * P        
//   - if Q != identity, then Q has order L (|E| = h2 * L, L prime)
//
// Output: accepted.csv
//   header = all columns from counted.csv + Bx_hex,By_hex:
//   id,p_hex,a_hex,d_hex,
//   order_hex,h2,L_hex,is_L_prime,bits_L,
//   twist_order_hex,twist_h2,twist_L_hex,twist_is_L_prime,twist_bits_L,
//   method,time_ms,
//   Bx_hex,By_hex

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "finite_field.h"
#include "edwards.h"
#include "ed_utils.h"

#ifndef BITS_MIN
#define BITS_MIN 1
#endif

#ifndef MAX_TRIES_GEN
#define MAX_TRIES_GEN 20000
#endif

#ifndef MAX_TRIES_POINT
#define MAX_TRIES_POINT 5000
#endif

static int str_eq(const char* a, const char* b){
    return strcmp(a,b)==0;
}

static int parse_u64(const char* s, uint64_t* out){
    if (!s) return 0;
    char* end = NULL;
    unsigned long long v = strtoull(s, &end, 10);
    if (s == end) return 0;
    *out = (uint64_t)v;
    return 1;
}

static void ff_set0(ff_elem* r){ memset(r->v, 0, sizeof r->v); }
static void ff_set1(ff_elem* r){ memset(r->v, 0, sizeof r->v); r->v[0] = 1; }
static void ff_set2(ff_elem* r){ memset(r->v, 0, sizeof r->v); r->v[0] = 2; }


static int is_identity(const ed_point_t* P, const ed_curve_t* EC){
    ff_elem x, y;
    if (!ed_to_affine(&x, &y, P, EC)) return 0;

    ff_elem zero, one;
    ff_set0(&zero);
    ff_set1(&one);

    return (memcmp(x.v, zero.v, sizeof x.v) == 0) &&
           (memcmp(y.v, one.v,  sizeof y.v) == 0);
}


// sample y; solve x^2 = (1 - y^2) / (a - d y^2)
static int sample_point_from_y(ed_point_t* P,
                               ff_elem* x_out,
                               ff_elem* y_out,
                               const ed_curve_t* EC,
                               const ff_ctx_t* F,
                               ed_rng_t* R)
{
    ff_elem y, y2, num, den, den_inv, x2, x;
    ff_elem one, zero;
    ff_set1(&one);
    ff_set0(&zero);

    for (int tries = 0; tries < MAX_TRIES_POINT; ++tries){
        // random y in Fp 
        ff_random(&y, R, F);
        ff_mul_e(&y2, &y, &y, F); // y^2

        // den = a - d*y^2 
        ff_elem d_y2;
        ff_mul_e(&d_y2, &EC->d, &y2, F);
        ff_sub_e(&den, &EC->a, &d_y2, F);

        // skip if den == 0 
        if (memcmp(den.v, zero.v, sizeof den.v) == 0)
            continue;

        // num = 1 - y^2 
        ff_sub_e(&num, &one, &y2, F);

        // x^2 = num / den 
        if (!ff_inv_e(&den_inv, &den, F)) continue;
        ff_mul_e(&x2, &num, &den_inv, F);

        // check if x2 is square in Fp
        if (!ff_is_square(&x2, F)) continue;
        if (!ff_sqrt_ts(&x, &x2, F)) continue;

        // convert (x,y) to extended coords
        if (!ed_from_affine(P, &x, &y, EC)) continue;

        if (x_out) *x_out = x;
        if (y_out) *y_out = y;
        return 1;
    }
    return 0;
}

int main(int argc, char** argv){
    const char* in_csv  = "data/counted.csv";
    const char* out_csv = "data/accepted.csv";
    int bits_min        = BITS_MIN;
    int verbose         = 0;

    for (int i = 1; i < argc; i++){
        if (str_eq(argv[i], "--in") && i+1 < argc){
            in_csv = argv[++i]; continue;
        }
        if (str_eq(argv[i], "--out") && i+1 < argc){
            out_csv = argv[++i]; continue;
        }
        if (str_eq(argv[i], "--bits-min") && i+1 < argc){
            bits_min = atoi(argv[++i]); continue;
        }
        if (str_eq(argv[i], "--verbose") && i+1 < argc){
            verbose = atoi(argv[++i]); continue;
        }
        fprintf(stderr,
            "Usage: %s [--in COUNTED.csv] [--out ACCEPTED.csv] "
            "[--bits-min 188] [--verbose 0|1]\n", argv[0]);
        return 1;
    }

    FILE* fin = fopen(in_csv, "r");
    if (!fin){
        perror("fopen input");
        return 1;
    }
    FILE* fout = fopen(out_csv, "w");
    if (!fout){
        perror("fopen output");
        fclose(fin);
        return 1;
    }

    // full header = counted.csv header + Bx_hex,By_hex
    fprintf(fout,
        "id,p_hex,a_hex,d_hex,"
        "order_hex,h2,L_hex,is_L_prime,bits_L,"
        "twist_order_hex,twist_h2,twist_L_hex,twist_is_L_prime,twist_bits_L,"
        "method,time_ms,"
        "Bx_hex,By_hex\n");

    char line[4096];
    if (!fgets(line, sizeof line, fin)){
        fclose(fin);
        fclose(fout);
        return 0;
    }

    int total    = 0;
    int accepted = 0;
    int row_idx  = 0;

    while (fgets(line, sizeof line, fin)){
        total++;
        row_idx++;

        char* cols[32];
        int   nc = 0;

        char* tok = strtok(line, ",\n\r");
        while (tok && nc < 32){
            cols[nc++] = tok;
            tok = strtok(NULL, ",\n\r");
        }

        if (nc < 16){
            if (verbose){
                fprintf(stderr, "[row %d] too few columns: %d\n", row_idx, nc);
            }
            continue;
        }

        // counted.csv header:
        const char* id                  = cols[0];
        const char* p_hex               = cols[1];
        const char* a_hex               = cols[2];
        const char* d_hex               = cols[3];
        const char* order_hex           = cols[4];
        const char* h2_str              = cols[5];
        const char* L_hex               = cols[6];
        const char* isL_str             = cols[7];
        const char* bitsL_str           = cols[8];
        const char* twist_order_hex     = cols[9];
        const char* twist_h2_str        = cols[10];
        const char* twist_L_hex         = cols[11];
        const char* twist_isL_str       = cols[12];
        const char* twist_bitsL_str     = cols[13];
        const char* method_str          = cols[14];
        const char* time_ms_str         = cols[15];

        int isL    = atoi(isL_str);
        int bits_L = atoi(bitsL_str);

        if (!isL){
            if (verbose) fprintf(stderr,"[id=%s] is_L_prime=0\n", id);
            continue;
        }
        if (bits_L < bits_min){
            if (verbose) fprintf(stderr,
                                 "[id=%s] bits_L=%d < %d\n",
                                 id, bits_L, bits_min);
            continue;
        }

        uint64_t h2;
        if (!parse_u64(h2_str, &h2)){
            if (verbose) fprintf(stderr,"[id=%s] bad h2\n", id);
            continue;
        }
        if (!(h2 == 1 || h2 == 4 || h2 == 8)){
            if (verbose) fprintf(stderr,
                                 "[id=%s] h2=%llu not in {1,4,8}\n",
                                 id, (unsigned long long)h2);
            continue;
        }

        ff_ctx_t F;
        ff_init_hex(&F, p_hex);

        ff_elem A, D, L;
        if (!ff_from_hex(a_hex, &A)){
            if (verbose) fprintf(stderr,"[id=%s] bad a_hex\n", id);
            continue;
        }
        if (!ff_from_hex(d_hex, &D)){
            if (verbose) fprintf(stderr,"[id=%s] bad d_hex\n", id);
            continue;
        }
        if (!ff_from_hex(L_hex, &L)){
            if (verbose) fprintf(stderr,"[id=%s] bad L_hex\n", id);
            continue;
        }

        ed_curve_t EC;
        if (!ed_init(&EC, &F, &A, &D)){
            if (verbose) fprintf(stderr,"[id=%s] ed_init failed\n", id);
            continue;
        }

        // RNG seeded
        ed_rng_t R;
        ed_rng_seed(&R, 0xC0FFEE123456789ULL ^ F.p[0] ^ (uint64_t)row_idx);

        ed_point_t P, Q;
        ff_elem x,y;
        int found = 0;

        for (int tries = 0; tries < MAX_TRIES_GEN && !found; ++tries){
            if (!sample_point_from_y(&P, &x, &y, &EC, &F, &R))
                continue;

            // cofactor clearing: Q = h2 * P
            if (h2 == 1){
                Q = P;
            } else {
                ed_scalar_mul_u64(&Q, &P, h2, &EC);
            }

            if (is_identity(&Q, &EC)){
                continue;
            }

            ff_elem Bx, By;
            ed_to_affine(&Bx, &By, &Q, &EC);

            char xb[2 + 16*FF_NLIMBS + 1];
            char yb[2 + 16*FF_NLIMBS + 1];
            ff_to_hex(&Bx, xb, sizeof xb);
            ff_to_hex(&By, yb, sizeof yb);

            fprintf(fout,
                "%s,%s,%s,%s,"
                "%s,%s,%s,%s,%s,"
                "%s,%s,%s,%s,%s,"
                "%s,%s,"
                "%s,%s\n",
                id, p_hex, a_hex, d_hex,
                order_hex, h2_str, L_hex, isL_str, bitsL_str,
                twist_order_hex, twist_h2_str, twist_L_hex, twist_isL_str, twist_bitsL_str,
                method_str, time_ms_str,
                xb, yb);

            found   = 1;
            accepted++;
            break;
        }

        if (!found && verbose){
            fprintf(stderr,"[id=%s] no generator found after %d tries\n",
                    id, MAX_TRIES_GEN);
        }
    }

    fclose(fin);
    fclose(fout);

    fprintf(stdout,
            "Wrote %d / %d accepted curves to %s\n",
            accepted, total, out_csv);
    return 0;
}
