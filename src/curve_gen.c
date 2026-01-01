#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>

#include "finite_field.h"
#include "edwards.h"
#include "ed_utils.h"
#include "cost.h"

#define MAX_TRIES 200000


static int pick_d(ff_elem* d, const ff_elem* a, const ff_ctx_t* F, ed_rng_t* R, int require_nonsquare) {
    // uniform sampling in Fp, with control over square/non-square distribution

    // choose a fixed non-residue v
    ff_elem nu;
    if (require_nonsquare) {
        for (;;) {
            ff_random(&nu, R, F);
            if (!ff_eq_zero(&nu) && ff_legendre(&nu, F) == -1)
                break;
        }
    }

    ff_elem t, dd;
    int tries = 0;
    for (;;) {
        if (++tries > MAX_TRIES)
            return 0;

        // Sample t from F_p uniformly
        ff_random(&t, R, F);
        if (ff_eq_zero(&t))
            continue;

        /* Construct d:
           - if require_nonsquare: d = v * t^2  (uniform over non-squares)
           - else:                 d = t^2      (uniform over squares)
        */
        ff_mul_e(&dd, &t, &t, F); // dd = t^2 
        if (require_nonsquare)
            ff_mul_e(&dd, &dd, &nu, F);

        // reject d == 0 or d == a
        if (ff_eq_zero(&dd))
            continue;
        int is_a = 1;
        for (int i = 0; i < FF_NLIMBS; i++) {
            if (dd.v[i] != a->v[i]) {
                is_a = 0;
                break;
            }
        }
        if (is_a)
            continue;

        *d = dd;
        return 1;
    }
}

static void write_csv_header(FILE* out) {
    fprintf(out, "p_hex,a_hex,d_hex,a_is_square,d_is_square,complete_ok\n");
}

static void write_csv_line(FILE* out, const ff_ctx_t* F, const ff_elem* a, const ff_elem* d, int complete_ok) {
    char p_hex[1024], a_hex[1024], d_hex[1024];

    // convert modulus p (limb array)
    ff_elem p_elem; memcpy(p_elem.v, F->p, sizeof p_elem.v);

    ff_to_hex(&p_elem, p_hex, sizeof(p_hex));
    ff_to_hex(a, a_hex, sizeof(a_hex));
    ff_to_hex(d, d_hex, sizeof(d_hex));

    int a_is_square = ff_legendre(a, F) == 1;
    int d_is_square = ff_legendre(d, F) == 1;

    fprintf(out, "%s,%s,%s,%d,%d,%d\n",
            p_hex, a_hex, d_hex, a_is_square, d_is_square, complete_ok);
}


int main(int argc, char** argv) {
    uint64_t count = 100;
    char p_hex[1024] = {0};
    char a_hex[1024] = {0};
    char out_path[1024] = "candidates.csv";
    int complete = 1;

    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--count") && i + 1 < argc) {
            count = strtoull(argv[++i], NULL, 10);
        } else if (!strcmp(argv[i], "--p") && i + 1 < argc) {
            strncpy(p_hex, argv[++i], sizeof(p_hex));
        } else if (!strcmp(argv[i], "--a") && i + 1 < argc) {
            strncpy(a_hex, argv[++i], sizeof(a_hex));
        } else if (!strcmp(argv[i], "--out") && i + 1 < argc) {
            strncpy(out_path, argv[++i], sizeof(out_path));
        } else if (!strcmp(argv[i], "--complete") && i + 1 < argc) {
            complete = atoi(argv[++i]);
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
        }
    }

    if (p_hex[0] == 0) {
        fprintf(stderr, "Usage: %s --p <hex> --a <hex> [--count N] [--complete 1] [--out file.csv]\n", argv[0]);
        return 1;
    }

    ff_ctx_t F;
    ff_init_hex(&F, p_hex);

    ff_elem A;
    if (a_hex[0]) {
        if (!ff_from_hex(a_hex, &A)) {
            fprintf(stderr, "Invalid a hex\n");
            return 1;
        }
        ff_red_e(&A, &F);
    } else {
        ff_from_u64(1, &A);
    }

    ed_rng_t R;
    ed_rng_seed(&R, (uint64_t)time(NULL) ^ 0xABCDEF123456789ULL);

    FILE* out = fopen(out_path, "w");
    if (!out) {
        perror("fopen");
        return 1;
    }
    write_csv_header(out);

    uint64_t ok = 0;
    for (uint64_t i = 0; i < count; i++) {
        ff_elem D;
        if (!pick_d(&D, &A, &F, &R, complete))
            continue;

        // completeness: a!=0, d!=0, a!=d, a is square, d is non-square 
        int comp_ok = ed_complete_candidate(&A, &D, &F);

        write_csv_line(out, &F, &A, &D, comp_ok);
        ok++;
    }

    fclose(out);
    printf("Generated %" PRIu64 " candidates written to %s\n", ok, out_path);
    return 0;
}
