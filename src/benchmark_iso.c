#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "finite_field.h"
#include "edwards.h"
#include "ed_utils.h"
#include "cost.h"

static const int CHAIN_KS[] = {10,20,50};
#define N_CHAIN_K (sizeof(CHAIN_KS)/sizeof(CHAIN_KS[0]))

static inline uint64_t ns_diff(const struct timespec* a,
                               const struct timespec* b)
{
    uint64_t s  = (uint64_t)(b->tv_sec - a->tv_sec);
    int64_t  ns = (int64_t)(b->tv_nsec - a->tv_nsec);
    return s*1000000000ull + (uint64_t)ns;
}

static volatile ff_elem sink;

static void iso4_construct(ff_elem* a2, ff_elem* d2,
                           const ff_elem* a,
                           const ff_elem* d,
                           const ff_ctx_t* F)
{
    COUNTED_NEG(a2, a, F);
    COUNTED_SUB(d2, d, a, F);

    sink = *a2;
    sink = *d2;
}

static void iso4_eval(ff_elem* x2, ff_elem* y2,
                      const ff_elem* x,
                      const ff_elem* y,
                      const ff_elem* a,
                      const ff_ctx_t* F)
{
    ff_elem xx, yy, ax2;
    ff_elem num, den;
    ff_elem two;

    ff_from_u64(2, &two);

    COUNTED_MUL(&xx, x, x, F);
    COUNTED_MUL(&yy, y, y, F);
    COUNTED_MUL(&ax2, a, &xx, F);

    COUNTED_MUL(&num, x, y, F);
    COUNTED_MUL(&num, &num, &two, F);
    COUNTED_SUB(&den, &yy, &ax2, F);
    ff_inv_e(&den, &den, F);
    COUNTED_MUL(x2, &num, &den, F);

    COUNTED_ADD(&num, &yy, &ax2, F);
    COUNTED_SUB(&den, &two, &yy, F);
    COUNTED_SUB(&den, &den, &ax2, F);
    ff_inv_e(&den, &den, F);
    COUNTED_MUL(y2, &num, &den, F);

    sink = *x2;
    sink = *y2;
}

static void iso4_kernel(ed_point_t* K,
                        const ed_curve_t* EC,
                        ed_rng_t* R)
{
    ed_point_t P, T2, T4;
    ff_elem x,y;

    if (!ed_sample_point_from_x(&P,&x,&y,EC,EC->f,R))
        return;

    ed_double(&T2, &P, EC);
    ed_double(&T4, &T2, EC);
    ed_double(K, &T4, EC);
}

static void bench_stage(FILE* fout,
                        const char* id,
                        const char* p_hex,
                        const char* a_hex,
                        const char* d_hex,
                        const char* stage,
                        int k,
                        int iters,
                        const ff_ctx_t* F,
                        const ed_curve_t* EC,
                        ed_rng_t* R)
{
    uint64_t t_ns[iters], Ms[iters], Ss[iters], As[iters];
    uint64_t sum_t=0,sumM=0,sumS=0,sumA=0;

    for (int it=0; it<iters; it++){

        ff_elem a = EC->a;
        ff_elem d = EC->d;

        ed_point_t P;
        ff_elem x,y;
        if (!ed_sample_point_from_x(&P,&x,&y,EC,F,R)){
            fprintf(stderr,"sample_point failed\n");
            return;
        }

        ff_cost_reset();
        struct timespec t0,t1;
        clock_gettime(CLOCK_MONOTONIC,&t0);

        if (!strcmp(stage,"kernel")){
            ed_point_t K;
            iso4_kernel(&K, EC, R);
        }
        else if (!strcmp(stage,"construct")){
            iso4_construct(&a,&d,&a,&d,F);
        }
        else if (!strcmp(stage,"eval")){
            iso4_eval(&x,&y,&x,&y,&a,F);
        }
        else if (!strcmp(stage,"chain")){
            for (int i=0;i<k;i++){
                iso4_construct(&a,&d,&a,&d,F);
                iso4_eval(&x,&y,&x,&y,&a,F);
            }
        }

        clock_gettime(CLOCK_MONOTONIC,&t1);
        ff_cost_t c = ff_cost_snapshot();

        uint64_t dt = ns_diff(&t0,&t1);
        t_ns[it]=dt; Ms[it]=c.M; Ss[it]=c.S; As[it]=c.A;
        sum_t+=dt; sumM+=c.M; sumS+=c.S; sumA+=c.A;
    }

    double avg_t=(double)sum_t/iters;
    double avgM=(double)sumM/iters;
    double avgS=(double)sumS/iters;
    double avgA=(double)sumA/iters;

    fprintf(fout,
        "%s,%s,%s,%s,%s,%d,%d",
        id,p_hex,a_hex,d_hex,stage,k,iters);

    for (int it=0; it<iters; it++){
        fprintf(fout,",%d,%llu,%llu,%llu,%llu",
            it+1,
            (unsigned long long)t_ns[it],
            (unsigned long long)Ms[it],
            (unsigned long long)Ss[it],
            (unsigned long long)As[it]);
    }

    fprintf(fout,
        ",%.3f,%.3f,%.3f,%.3f,%llu,%llu,%llu\n",
        avg_t,avgM,avgS,avgA,
        (unsigned long long)sumM,
        (unsigned long long)sumS,
        (unsigned long long)sumA);
}

int main(int argc, char** argv)
{
    const char* in="data/accepted.csv";
    const char* out="data/bench_iso.csv";
    int iters=10;

    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"--in") && i+1<argc) in=argv[++i];
        else if (!strcmp(argv[i],"--out") && i+1<argc) out=argv[++i];
        else if (!strcmp(argv[i],"--iters") && i+1<argc) iters=atoi(argv[++i]);
    }

    FILE* fin=fopen(in,"r");
    FILE* fout=fopen(out,"w");
    if (!fin||!fout){ perror("fopen"); return 1; }

    fprintf(fout,
        "id,p_hex,a_hex,d_hex,stage,k,N");

    for (int i=1;i<=iters;i++)
        fprintf(fout,",iter%d,time%d_ns,M%d,S%d,A%d",i,i,i,i,i);

    fprintf(fout,
        ",ns_avg,M_avg,S_avg,A_avg,M_total,S_total,A_total\n");

    char line[4096];
    fgets(line,sizeof line,fin);

    while (fgets(line,sizeof line,fin)){
        char* c[8]; int n=0;
        char* t=strtok(line,",\n\r");
        while(t&&n<8){ c[n++]=t; t=strtok(NULL,",\n\r"); }
        if (n<4) continue;

        const char* id=c[0];
        const char* p_hex=c[1];
        const char* a_hex=c[2];
        const char* d_hex=c[3];

        ff_ctx_t F; ff_init_hex(&F,p_hex);
        ff_elem a0,d0;
        ff_from_hex(a_hex,&a0);
        ff_from_hex(d_hex,&d0);

        ed_curve_t EC;
        ed_init(&EC,&F,&a0,&d0);

        ed_rng_t R;
        ed_rng_seed(&R, 0xC0FFEE ^ F.p[0]);

        bench_stage(fout,id,p_hex,a_hex,d_hex,
                    "kernel",1,iters,&F,&EC,&R);
        bench_stage(fout,id,p_hex,a_hex,d_hex,
                    "construct",1,iters,&F,&EC,&R);
        bench_stage(fout,id,p_hex,a_hex,d_hex,
                    "eval",1,iters,&F,&EC,&R);

        for (size_t i=0;i<N_CHAIN_K;i++)
            bench_stage(fout,id,p_hex,a_hex,d_hex,
                        "chain",CHAIN_KS[i],iters,&F,&EC,&R);
    }

    fclose(fin);
    fclose(fout);
    return 0;
}
