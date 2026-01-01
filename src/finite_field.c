#include "finite_field.h"
#include <string.h>
#include <assert.h>

static int ff_clamp_nlimbs(const uint64_t a[FF_NLIMBS]) {
    int n = FF_NLIMBS;
    while (n > 1 && a[n - 1] == 0) --n;
    return n;
}

static inline uint64_t add_carry64(uint64_t a, uint64_t b, uint64_t *carry) {
    __uint128_t s = (__uint128_t)a + b + *carry;
    *carry = (uint64_t)(s >> 64);
    return (uint64_t)s;
}

static inline uint64_t sub_borrow64(uint64_t a, uint64_t b, uint64_t *borrow) {
    __uint128_t d = (__uint128_t)a - b - *borrow;
    *borrow = (uint64_t)((d >> 127) & 1);
    return (uint64_t)d;
}

static int cmp_limbs(const uint64_t *a, const uint64_t *b, int n) {
    for (int i = n - 1; i >= 0; --i) {
        if (a[i] != b[i]) return (a[i] > b[i]) ? 1 : -1;
    }
    return 0;
}

static void cpy_limbs(uint64_t *dst, const uint64_t *src) {
    for (int i = 0; i < FF_NLIMBS; ++i) dst[i] = src[i];
}

static void clr_limbs(uint64_t *a) {
    for (int i = 0; i < FF_NLIMBS; ++i) a[i] = 0;
}


// hex pairing + format
static int hex_val(int c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

bool ff_from_hex(const char *hex, ff_elem *a) {
    if (!hex || !a) return false;
    clr_limbs(a->v);
    if (hex[0] == '0' && (hex[1] == 'x' || hex[1] == 'X')) hex += 2;
    size_t len = strlen(hex);
    if (!len) return false;

    size_t nib = 0;
    for (ssize_t i = (ssize_t)len - 1; i >= 0; --i, ++nib) {
        int v = hex_val(hex[i]);
        if (v < 0) return false;
        size_t bit  = nib * 4;
        size_t limb = bit / 64;
        size_t off  = bit % 64;
        if (limb >= (size_t)FF_NLIMBS) return false;
        a->v[limb] |= (uint64_t)v << off;
        if (off > 60 && limb + 1 < FF_NLIMBS)
            a->v[limb + 1] |= (uint64_t)v >> (64 - off);
    }
    return true;
}

static void hex_append64(char *dst, size_t dst_sz, size_t *pos, uint64_t limb, int keep_all) {
    static const char *H = "0123456789abcdef";
    char buf[16];
    for (int i = 15; i >= 0; --i) { buf[i] = H[limb & 0xF]; limb >>= 4; }
    int start = 0;
    if (!keep_all) {
        while (start < 16 && buf[start] == '0') ++start;
        if (start == 16) { buf[15] = '0'; start = 15; }
    }
    for (int i = start; i < 16 && *pos + 1 < dst_sz; ++i) dst[(*pos)++] = buf[i];
}

void ff_to_hex(const ff_elem *a, char *out, size_t sz) {
    if (!out || sz < 3) return;
    size_t pos = 0; out[pos++] = '0'; out[pos++] = 'x';
    int seen = 0;
    for (int i = FF_NLIMBS - 1; i >= 0; --i) {
        if (!seen && a->v[i] == 0 && i > 0) continue;
        hex_append64(out, sz, &pos, a->v[i], seen ? 1 : 0);
        seen = 1;
    }
    if (!seen) out[pos++] = '0';
    out[pos] = 0;
}

void ff_from_u64(uint64_t x, ff_elem *a) {
    clr_limbs(a->v);
    a->v[0] = x;
}

bool ff_to_u64(const ff_elem *a, uint64_t *out) {
    for (int i = 1; i < FF_NLIMBS; ++i) if (a->v[i]) return false;
    if (out) *out = a->v[0];
    return true;
}

static uint64_t inv64_mod2_64(uint64_t x) {
    uint64_t inv = 1;
    for (int i = 0; i < 6; ++i) inv *= 2 - x * inv;
    return inv;
}

void ff_red_e(ff_elem *a, const ff_ctx_t *ctx);


static void compute_R2(uint64_t R2[FF_NLIMBS], const ff_ctx_t *ctx) {
    ff_elem x; clr_limbs(x.v); x.v[0] = 1;
    int bits = FF_NLIMBS * 64 * 2; // R^2 = 2^(128*FF_NLIMBS)
    for (int i = 0; i < bits; ++i) {
        uint64_t carry = 0, next;
        for (int j = 0; j < FF_NLIMBS; ++j) {
            next = x.v[j] >> 63;
            x.v[j] = (x.v[j] << 1) | carry;
            carry = next;
        }
        ff_red_e(&x, ctx);
    }
    cpy_limbs(R2, x.v);
}

bool ff_init(ff_ctx_t *ctx, uint64_t p_u64) {
    if (!ctx) return false;
    if (p_u64 < 3ull || !(p_u64 & 1ull)) return false;
    clr_limbs(ctx->p); ctx->p[0] = p_u64;
    ctx->nlimbs = ff_clamp_nlimbs(ctx->p);
    uint64_t inv = inv64_mod2_64(ctx->p[0]);
    ctx->n0 = (uint64_t)(0 - inv);
    compute_R2(ctx->R2, ctx);
    return true;
}

bool ff_init_hex(ff_ctx_t *ctx, const char *hex_p) {
    if (!ctx) return false;
    ff_elem pe;
    if (!ff_from_hex(hex_p, &pe)) return false;
    if (!(pe.v[0] & 1ull)) return false;
    int n = ff_clamp_nlimbs(pe.v);
    if (n == 1 && pe.v[0] < 3ull) return false;

    cpy_limbs(ctx->p, pe.v);
    ctx->nlimbs = ff_clamp_nlimbs(ctx->p);
    uint64_t inv = inv64_mod2_64(ctx->p[0]);
    ctx->n0 = (uint64_t)(0 - inv);
    compute_R2(ctx->R2, ctx);
    return true;
}

// Montgomery MUL(CIOS), conversion
static void mont_mul(uint64_t r[FF_NLIMBS], const uint64_t a[FF_NLIMBS],
                     const uint64_t b[FF_NLIMBS], const ff_ctx_t *ctx) {
    uint64_t t[FF_NLIMBS + 1];
    for (int i = 0; i <= FF_NLIMBS; ++i) t[i] = 0;

    for (int i = 0; i < FF_NLIMBS; ++i) {
        __uint128_t carry = 0;
        for (int j = 0; j < FF_NLIMBS; ++j) {
            __uint128_t u = (__uint128_t)a[j] * b[i] + t[j] + carry;
            t[j] = (uint64_t)u;
            carry = u >> 64;
        }
        t[FF_NLIMBS] = (uint64_t)((__uint128_t)t[FF_NLIMBS] + carry);

        uint64_t m = t[0] * ctx->n0;

        __uint128_t c2 = 0, sum;
        for (int j = 0; j < FF_NLIMBS; ++j) {
            sum = (__uint128_t)t[j] + (__uint128_t)m * ctx->p[j] + c2;
            if (j > 0) t[j - 1] = (uint64_t)sum;
            c2 = sum >> 64;
        }
        __uint128_t top = (__uint128_t)t[FF_NLIMBS] + c2;
        t[FF_NLIMBS - 1] = (uint64_t)top;
        t[FF_NLIMBS]     = (uint64_t)(top >> 64);
    }
    for (int i = 0; i < FF_NLIMBS; ++i) r[i] = t[i];

    uint64_t tmp[FF_NLIMBS], br = 0;
    for (int i = 0; i < FF_NLIMBS; ++i) tmp[i] = sub_borrow64(r[i], ctx->p[i], &br);
    uint64_t mask = (uint64_t)-(uint64_t)(br == 0);
    for (int i = 0; i < FF_NLIMBS; ++i) r[i] = (tmp[i] & mask) | (r[i] & ~mask);
}

static void to_mont(uint64_t r[FF_NLIMBS], const uint64_t a[FF_NLIMBS], const ff_ctx_t *ctx) {
    mont_mul(r, a, ctx->R2, ctx);
}

static void from_mont(uint64_t r[FF_NLIMBS], const uint64_t aM[FF_NLIMBS], const ff_ctx_t *ctx) {
    uint64_t one[FF_NLIMBS] = {0}; one[0] = 1;
    mont_mul(r, aM, one, ctx);
}

// 512 arthmetic API
void ff_red_e(ff_elem *a, const ff_ctx_t *ctx) {
    uint64_t tmp[FF_NLIMBS];
    uint64_t br = 0;
    for (int i = 0; i < FF_NLIMBS; ++i) tmp[i] = sub_borrow64(a->v[i], ctx->p[i], &br);
    uint64_t mask = (uint64_t)-(uint64_t)(br == 0);
    for (int i = 0; i < FF_NLIMBS; ++i) a->v[i] = (tmp[i] & mask) | (a->v[i] & ~mask);

    br = 0;
    for (int i = 0; i < FF_NLIMBS; ++i) tmp[i] = sub_borrow64(a->v[i], ctx->p[i], &br);
    mask = (uint64_t)-(uint64_t)(br == 0);
    for (int i = 0; i < FF_NLIMBS; ++i) a->v[i] = (tmp[i] & mask) | (a->v[i] & ~mask);
}

void ff_add_e(ff_elem *r, const ff_elem *a, const ff_elem *b, const ff_ctx_t *ctx) {
    uint64_t c = 0;
    for (int i = 0; i < FF_NLIMBS; ++i) r->v[i] = add_carry64(a->v[i], b->v[i], &c);
    ff_red_e(r, ctx);
}

void ff_sub_e(ff_elem *r, const ff_elem *a, const ff_elem *b, const ff_ctx_t *ctx) {
    uint64_t br = 0;
    for (int i = 0; i < FF_NLIMBS; ++i) r->v[i] = sub_borrow64(a->v[i], b->v[i], &br);
    if (br) {
        uint64_t c = 0;
        for (int i = 0; i < FF_NLIMBS; ++i) r->v[i] = add_carry64(r->v[i], ctx->p[i], &c);
    }
}

void ff_neg_e(ff_elem *r, const ff_elem *a, const ff_ctx_t *ctx) {
    uint64_t zero[FF_NLIMBS] = {0};
    if (cmp_limbs(a->v, zero, FF_NLIMBS) == 0) { clr_limbs(r->v); return; }
    uint64_t br = 0;
    for (int i = 0; i < FF_NLIMBS; ++i) r->v[i] = sub_borrow64(ctx->p[i], a->v[i], &br);
}

void ff_mul_e(ff_elem *r, const ff_elem *a, const ff_elem *b, const ff_ctx_t *ctx) {
    uint64_t aa[FF_NLIMBS], bb[FF_NLIMBS], aM[FF_NLIMBS], bM[FF_NLIMBS], cM[FF_NLIMBS];
    cpy_limbs(aa, a->v); ff_red_e((ff_elem *)aa, ctx);
    cpy_limbs(bb, b->v); ff_red_e((ff_elem *)bb, ctx);
    to_mont(aM, aa, ctx);
    to_mont(bM, bb, ctx);
    mont_mul(cM, aM, bM, ctx);
    from_mont(r->v, cM, ctx);
}

void ff_pow_e(ff_elem *r, const ff_elem *a, const ff_elem *e, const ff_ctx_t *ctx) {
    uint64_t baseM[FF_NLIMBS], resM[FF_NLIMBS];
    uint64_t one[FF_NLIMBS] = {1};
    uint64_t aa[FF_NLIMBS]; cpy_limbs(aa, a->v); ff_red_e((ff_elem *)aa, ctx);
    to_mont(baseM, aa, ctx);
    to_mont(resM, one, ctx);

    for (int i = 0; i < FF_NLIMBS; ++i) {
        uint64_t w = e->v[i];
        for (int b = 0; b < 64; ++b) {
            if (w & 1ull) mont_mul(resM, resM, baseM, ctx);
            w >>= 1;
            mont_mul(baseM, baseM, baseM, ctx);
        }
    }
    from_mont(r->v, resM, ctx);
}

bool ff_inv_e(ff_elem *r, const ff_elem *a, const ff_ctx_t *ctx) {
    uint64_t zero[FF_NLIMBS] = {0};
    if (cmp_limbs(a->v, zero, FF_NLIMBS) == 0) { clr_limbs(r->v); return false; }
    // Fermat's little theorem for odd p
    ff_elem exp; cpy_limbs(exp.v, ctx->p);
    uint64_t br = 0;
    exp.v[0] = sub_borrow64(exp.v[0], 2, &br);
    for (int i = 1; i < FF_NLIMBS; ++i) exp.v[i] = sub_borrow64(exp.v[i], 0, &br);
    ff_pow_e(r, a, &exp, ctx);
    return true;
}

// 64bit API
static void ensure_u64_ctx(const ff_ctx_t *ctx) {
    assert(ctx->nlimbs == 1 && "64-bit API used with modulus > 64 bits");
}

uint64_t ff_red(uint64_t a, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    return a % ctx->p[0];
}

uint64_t ff_add(uint64_t a, uint64_t b, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    uint64_t p = ctx->p[0];
    uint64_t s = a + b;
    uint64_t mask = (uint64_t)-((__int128)s - (__int128)p >= 0);
    return s - (p & mask);
}

uint64_t ff_sub(uint64_t a, uint64_t b, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    uint64_t p = ctx->p[0];
    uint64_t d = a - b;
    uint64_t mask = (uint64_t)-((a < b) & 1u);
    return d + (p & mask);
}

uint64_t ff_neg(uint64_t a, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    uint64_t p = ctx->p[0];
    uint64_t mask = (uint64_t)-(uint64_t)(a != 0);
    return (p - a) & mask;
}

uint64_t ff_mul(uint64_t a, uint64_t b, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    ff_elem ae, be, re;
    ff_from_u64(a, &ae);
    ff_from_u64(b, &be);
    ff_mul_e(&re, &ae, &be, ctx);
    uint64_t out = 0; (void)ff_to_u64(&re, &out);
    return out;
}

uint64_t ff_pow(uint64_t a, uint64_t e, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    ff_elem ae, ee, re;
    ff_from_u64(a, &ae);
    ff_from_u64(0, &ee); ee.v[0] = e;
    ff_pow_e(&re, &ae, &ee, ctx);
    uint64_t out = 0; (void)ff_to_u64(&re, &out);
    return out;
}

uint64_t ff_inv(uint64_t a, const ff_ctx_t *ctx) {
    ensure_u64_ctx(ctx);
    if (a == 0) return 0;
    ff_elem ae, re;
    ff_from_u64(a, &ae);
    if (!ff_inv_e(&re, &ae, ctx)) return 0;
    uint64_t out = 0; (void)ff_to_u64(&re, &out);
    return out;
}
