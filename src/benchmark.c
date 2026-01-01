// benchmark.c
//
// Benchmarks different SMUL (binary, w-NAF, fixed-window, ladder) / ADD / DOUBLE operations for curves from accepted.csv
//
// Input: accepted.csv
// Header (accepted.csv):
//    id,p_hex,a_hex,d_hex,order_hex,h2,L_hex,is_L_prime,bits_L,
//    twist_order_hex,twist_h2,twist_L_hex,twist_is_L_prime,twist_bits_L,
//    method,time_ms,Bx_hex,By_hex
//
// Output (bench.csv):
//   id,p_hex,a_hex,d_hex,h,L_hex,Bx_hex,By_hex,
//   algo,win,nbits,scalar_hex,N,
//   iter1,time1_ns,M1,S1,A1,
//   iter2,time2_ns,M2,S2,A2,
//   ...
//   iterN,timeN_ns,MN,SN,AN,
//   ns_avg,M_avg,S_avg,A_avg,M_total,S_total,A_total

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "finite_field.h"
#include "edwards.h"
#include "ed_utils.h"
#include "cost.h"
#include "edwards_acc.h"  


#ifndef DEFAULT_ITERS
#define DEFAULT_ITERS 20
#endif

// windows size
static const int WNAF_WS[]  = {3, 4, 5, 6,7};
static const size_t N_WNAF_W = sizeof(WNAF_WS)/sizeof(WNAF_WS[0]);

static const int FIXED_WS[] = {3,4, 5, 6,7};
static const size_t N_FIXED_W = sizeof(FIXED_WS)/sizeof(FIXED_WS[0]);

static inline uint64_t ns_diff(const struct timespec* t0,
                               const struct timespec* t1)
{
    uint64_t s  = (uint64_t)(t1->tv_sec  - t0->tv_sec);
    int64_t  ns = (int64_t)(t1->tv_nsec - t0->tv_nsec);
    return s*1000000000ull + (uint64_t)ns;
}

static int parse_u64_dec(const char* s, uint64_t* out){
    if (!s) return 0;
    char* end = NULL;
    unsigned long long v = strtoull(s,&end,10);
    if (s == end) return 0;
    *out = (uint64_t)v;
    return 1;
}

static size_t ff_to_bits_lsb_first(const ff_elem *x,
                                   uint8_t *bits,
                                   size_t bits_cap)
{
    const size_t total_bits = (size_t)FF_NLIMBS * 64u;
    if (bits_cap < total_bits) return 0;

    memset(bits, 0, (total_bits + 7) / 8);

    for (size_t limb = 0; limb < (size_t)FF_NLIMBS; ++limb) {
        uint64_t w = x->v[limb];
        for (int b = 0; b < 64; ++b){
            if ((w >> b) & 1u){
                size_t bit_index = limb*64u + (size_t)b;
                bits[bit_index>>3] |= (uint8_t)(1u << (bit_index & 7));
            }
        }
    }

    size_t nbits = total_bits;
    while (nbits > 0){
        size_t i = nbits - 1;
        uint8_t byte = bits[i>>3];
        if (((byte >> (i & 7)) & 1u) == 0) nbits--;
        else break;
    }
    return nbits;
}


typedef struct {
    ed_curve_t EC;
    ed_point_t B; // curve subgroup G generator (#G=L)
    ed_point_t P, Q, R;
    uint8_t*   k_bits;
    size_t     k_nbits;
    ed_fixed_precomp_t PC;  
    int wnaf_w;
    int fixed_w;
} bench_ctx_t;


static void ed_scalar_mul_binary_naive_bits(ed_point_t* R,
                                            const ed_point_t* P,
                                            const uint8_t* k_le_bits,
                                            size_t nbits,
                                            const ed_curve_t* ec)
{
    // R = identity (affine)
    ff_elem x0, y1;
    ff_from_u64(0,&x0);
    ff_from_u64(1,&y1);
    ed_from_affine(R, &x0, &y1, ec);

    for (size_t i = nbits; i-- > 0; ){
        ed_double(R, R, ec);
        uint8_t b = (k_le_bits[i>>3] >> (i & 7)) & 1u;
        if (b) ed_add(R, R, P, ec);
    }
}


static void once_add(void* vctx){
    bench_ctx_t* c = (bench_ctx_t*)vctx;
    ed_add(&c->R, &c->P, &c->Q, &c->EC);
    c->P = c->R;
}

static void once_double(void* vctx){
    bench_ctx_t* c = (bench_ctx_t*)vctx;
    ed_double(&c->R, &c->P, &c->EC);
    c->P = c->R;
}

// ladder 
static void once_smul_ladder(void* vctx){
    bench_ctx_t* c = (bench_ctx_t*)vctx;
    ed_scalar_mul_bits(&c->R, &c->B, c->k_bits, c->k_nbits, &c->EC);
}

// binary naive 
static void once_smul_binary_naive(void* vctx){
    bench_ctx_t* c = (bench_ctx_t*)vctx;
    ed_scalar_mul_binary_naive_bits(&c->R, &c->B, c->k_bits, c->k_nbits, &c->EC);
}

// w-NAF with window w in ctx->wnaf_w
static void once_smul_wnaf(void* vctx){
    bench_ctx_t* c = (bench_ctx_t*)vctx;
    ed_scalar_mul_wnaf_bits(&c->R, &c->B,
                            c->k_bits, c->k_nbits,
                            c->wnaf_w, &c->EC);
}

// fixed-windows with PC array
static void once_smul_fixed(void* vctx){
    bench_ctx_t* c = (bench_ctx_t*)vctx;
    ed_scalar_mul_fixed_bits(&c->R,
                             c->k_bits, c->k_nbits,
                             &c->PC, &c->EC);
}


static void bench_method_row(FILE* fout,
                             const char* id,
                             const char* p_hex,
                             const char* a_hex,
                             const char* d_hex,
                             const char* h_str,
                             const char* L_hex,
                             const char* Bx_hex,
                             const char* By_hex,
                             const char* algo,
                             int win,
                             size_t nbits_scalar,
                             const char* scalar_hex,
                             int N,
                             void (*once)(void*),
                             bench_ctx_t* ctx)
{
    if (N <= 0) return;

    uint64_t* t_ns = (uint64_t*)calloc((size_t)N, sizeof(uint64_t));
    uint64_t* Ms   = (uint64_t*)calloc((size_t)N, sizeof(uint64_t));
    uint64_t* Ss   = (uint64_t*)calloc((size_t)N, sizeof(uint64_t));
    uint64_t* As   = (uint64_t*)calloc((size_t)N, sizeof(uint64_t));
    if (!t_ns || !Ms || !Ss || !As){
        free(t_ns); free(Ms); free(Ss); free(As);
        fprintf(stderr,"[id=%s algo=%s] malloc failed – skipping\n", id, algo);
        return;
    }

    struct timespec t0,t1;
    uint64_t sum_ns = 0;
    uint64_t sumM = 0, sumS = 0, sumA = 0;

    for (int j=0; j<N; ++j){
        ff_cost_reset();
        ff_cost_t before = ff_cost_snapshot();

        clock_gettime(CLOCK_MONOTONIC, &t0);
        once(ctx);
        clock_gettime(CLOCK_MONOTONIC, &t1);

        ff_cost_t after = ff_cost_snapshot();

        uint64_t dt = ns_diff(&t0,&t1);
        uint64_t dM = after.M;
        uint64_t dS = after.S;
        uint64_t dA = after.A;

        t_ns[j] = dt;
        Ms[j]   = dM;
        Ss[j]   = dS;
        As[j]   = dA;

        sum_ns += dt;
        sumM   += dM;
        sumS   += dS;
        sumA   += dA;
    }

    double ns_avg = (double)sum_ns / (double)N;
    double M_avg  = (double)sumM  / (double)N;
    double S_avg  = (double)sumS  / (double)N;
    double A_avg  = (double)sumA  / (double)N;


    fprintf(fout,
        "%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%zu,%s,%d",
        id, p_hex, a_hex, d_hex, h_str, L_hex, Bx_hex, By_hex,
        algo, win, nbits_scalar, scalar_hex, N);

    for (int j=0; j<N; ++j){
        int iter_index = j+1;
        fprintf(fout,
                ",%d,%llu,%llu,%llu,%llu",
                iter_index,
                (unsigned long long)t_ns[j],
                (unsigned long long)Ms[j],
                (unsigned long long)Ss[j],
                (unsigned long long)As[j]);
    }

    fprintf(fout,
            ",%.6f,%.4f,%.4f,%.4f,%llu,%llu,%llu\n",
            ns_avg, M_avg, S_avg, A_avg,
            (unsigned long long)sumM,
            (unsigned long long)sumS,
            (unsigned long long)sumA);

    free(t_ns);
    free(Ms);
    free(Ss);
    free(As);
}


int main(int argc, char** argv){
    const char* in_path  = "data/accepted.csv";
    const char* out_path = "data/bench.csv";
    int iters = DEFAULT_ITERS;

    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"--in") && i+1<argc){
            in_path = argv[++i];
            continue;
        }
        if (!strcmp(argv[i],"--out") && i+1<argc){
            out_path = argv[++i];
            continue;
        }
        if (!strcmp(argv[i],"--iters") && i+1<argc){
            iters = atoi(argv[++i]);
            continue;
        }
        fprintf(stderr,
            "Usage: %s [--in ACCEPTED.csv] [--out BENCH.csv] [--iters N]\n",
            argv[0]);
        return 1;
    }

    if (iters <= 0){
        fprintf(stderr,"ERROR: iters must be > 0\n");
        return 1;
    }

    FILE* fin = fopen(in_path, "r");
    if (!fin){
        perror("fopen input");
        return 1;
    }
    FILE* fout = fopen(out_path, "w");
    if (!fout){
        perror("fopen output");
        fclose(fin);
        return 1;
    }

    fprintf(fout,
        "id,p_hex,a_hex,d_hex,h,L_hex,Bx_hex,By_hex,"
        "algo,win,nbits,scalar_hex,N");

    for (int j=1;j<=iters;j++){
        fprintf(fout,",iter%d,time%d_ns,M%d,S%d,A%d", j,j,j,j,j);
    }
    fprintf(fout,
        ",ns_avg,M_avg,S_avg,A_avg,M_total,S_total,A_total\n");

    char line[4096];
    if (!fgets(line, sizeof line, fin)){
        fclose(fin);
        fclose(fout);
        return 0;
    }

    int row_idx = 0;
    int bench_curves = 0;

    while (fgets(line, sizeof line, fin)){
        row_idx++;

        char* cols[32];
        int   nc = 0;
        char* tok = strtok(line, ",\n\r");
        while (tok && nc<32){
            cols[nc++] = tok;
            tok = strtok(NULL, ",\n\r");
        }

        if (nc < 18){
            fprintf(stderr,"[row %d] too few columns: %d (need >=18)\n", row_idx, nc);
            continue;
        }

        const char* id          = cols[0];
        const char* p_hex       = cols[1];
        const char* a_hex       = cols[2];
        const char* d_hex       = cols[3];
        const char* h_str       = cols[5];   
        const char* L_hex       = cols[6];
        const char* isL_str     = cols[7];
        const char* bitsL_str   = cols[8];
        const char* Bx_hex      = cols[16];
        const char* By_hex      = cols[17];

        int isL = atoi(isL_str);
        if (!isL){
            fprintf(stderr,"[id=%s] is_L_prime=0 – skip\n", id);
            continue;
        }

        int bitsL = atoi(bitsL_str);
        if (bitsL <= 0){
            fprintf(stderr,"[id=%s] bits_L<=0 – skip\n", id);
            continue;
        }

        uint64_t hval = 0;
        if (!parse_u64_dec(h_str, &hval)){
            fprintf(stderr,"[id=%s] bad h='%s' – skip\n", id, h_str);
            continue;
        }

        ff_ctx_t F;
        ff_init_hex(&F, p_hex);

        ff_elem A, D, Bx, By, L;
        if (!ff_from_hex(a_hex, &A)){
            fprintf(stderr,"[id=%s] bad a_hex\n", id);
            continue;
        }
        if (!ff_from_hex(d_hex, &D)){
            fprintf(stderr,"[id=%s] bad d_hex\n", id);
            continue;
        }
        if (!ff_from_hex(Bx_hex, &Bx)){
            fprintf(stderr,"[id=%s] bad Bx_hex\n", id);
            continue;
        }
        if (!ff_from_hex(By_hex, &By)){
            fprintf(stderr,"[id=%s] bad By_hex\n", id);
            continue;
        }
        if (!ff_from_hex(L_hex, &L)){
            fprintf(stderr,"[id=%s] bad L_hex\n", id);
            continue;
        }

        ed_curve_t EC;
        if (!ed_init(&EC, &F, &A, &D)){
            fprintf(stderr,"[id=%s] ed_init failed – skip\n", id);
            continue;
        }

        if (!ed_oncurve_affine(&Bx, &By, &EC)){
            fprintf(stderr,"[id=%s] base point NOT on curve – skip\n", id);
            continue;
        }

        ed_point_t B;
        ed_from_affine(&B, &Bx, &By, &EC);

        size_t bits_cap  = (size_t)FF_NLIMBS * 64u;
        size_t bytes_cap = (bits_cap + 7)/8;
        uint8_t* k_bits = (uint8_t*)calloc(bytes_cap, 1);
        if (!k_bits){
            fprintf(stderr,"[id=%s] calloc k_bits failed\n", id);
            continue;
        }
        size_t k_nbits = ff_to_bits_lsb_first(&L, k_bits, bits_cap);
        if (k_nbits == 0){
            fprintf(stderr,"[id=%s] k_nbits=0 – skip\n", id);
            free(k_bits);
            continue;
        }

        const char* scalar_hex  = L_hex;
        size_t nbits_scalar     = k_nbits;

        bench_ctx_t ctx;
        memset(&ctx, 0, sizeof ctx);
        ctx.EC      = EC;
        ctx.B       = B;
        ctx.P       = B;
        ctx.Q       = B;
        ctx.R       = B;
        ctx.k_bits  = k_bits;
        ctx.k_nbits = k_nbits;
        ctx.PC.tbl  = NULL;
        ctx.wnaf_w  = 0;
        ctx.fixed_w = 0;

        // add(B,B)
        ctx.P = ctx.B;
        ctx.Q = ctx.B;
        bench_method_row(fout,
            id,p_hex,a_hex,d_hex,h_str,L_hex,Bx_hex,By_hex,
            "add",0,nbits_scalar,scalar_hex,
            iters, once_add,&ctx);

        // double(B)
        ctx.P = ctx.B;
        bench_method_row(fout,
            id,p_hex,a_hex,d_hex,h_str,L_hex,Bx_hex,By_hex,
            "double",0,nbits_scalar,scalar_hex,
            iters, once_double,&ctx);

        // ladder
        bench_method_row(fout,
            id,p_hex,a_hex,d_hex,h_str,L_hex,Bx_hex,By_hex,
            "ladder",0,nbits_scalar,scalar_hex,
            iters, once_smul_ladder,&ctx);

        // binary naive 
        bench_method_row(fout,
            id,p_hex,a_hex,d_hex,h_str,L_hex,Bx_hex,By_hex,
            "binary",0,nbits_scalar,scalar_hex,
            iters, once_smul_binary_naive,&ctx);

        // w-NAF(w,B)
        for (size_t iw=0;iw<N_WNAF_W;iw++){
            ctx.wnaf_w = WNAF_WS[iw];
            bench_method_row(fout,
                id,p_hex,a_hex,d_hex,h_str,L_hex,Bx_hex,By_hex,
                "wNAF",ctx.wnaf_w,nbits_scalar,scalar_hex,
                iters, once_smul_wnaf,&ctx);
        }

        // fixed-window(w,B)
        for (size_t iw=0;iw<N_FIXED_W;iw++){
            ctx.fixed_w = FIXED_WS[iw];

            ed_fixed_precomp_t pc;
            memset(&pc, 0, sizeof pc);
            if (!ed_fixed_precompute(&pc, &ctx.B, ctx.fixed_w, &ctx.EC)){
                fprintf(stderr,"[id=%s] fixed precompute failed for w=%d – skip fixed\n",
                        id, ctx.fixed_w);
                continue;
            }
            ctx.PC = pc;

            bench_method_row(fout,
                id,p_hex,a_hex,d_hex,h_str,L_hex,Bx_hex,By_hex,
                "fixed",ctx.fixed_w,nbits_scalar,scalar_hex,
                iters, once_smul_fixed,&ctx);

            ed_fixed_precomp_free(&ctx.PC);
            memset(&ctx.PC,0,sizeof ctx.PC);
        }

        free(k_bits);
        bench_curves++;
    }

    fclose(fin);
    fclose(fout);

    fprintf(stderr, "Benchmark done for %d curves -> %s\n",
            bench_curves, out_path);

    return 0;
}
