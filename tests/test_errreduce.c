/* test_errreduce.c — cross-checks for the errreduce module (bead aic-t81):
 * cor_improvement / error reduction (approximate_algebras.tex:1317-1319). A THIN
 * wrapper over dhom (§8), so the tests focus on the cor_improvement CONTRACT:
 * a delta-inclusion in, a CERTIFIED c_0*eps-inclusion out, with a DIMENSION-
 * INDEPENDENT measured c_0 (the .tex:484/461 universality canary, the headline).
 *
 * Cross-check ladder (CLAUDE.md Rule 6/7; mutation-proven, per-test comments):
 *   T1 eta=0 ORACLE — A=M_n identity channel, v(X)=U X U^dag an exact *-iso (an
 *      exact 0-inclusion). errreduce returns vtilde=v unchanged, ALL four
 *      inclusion-defects machine-zero, c_0*eps = 0.
 *      MUTATION: skip lem_approx and inject a non-inclusion (collapse a basis op)
 *      -> the input-delta-inclusion guard aborts (verified RED, see comment).
 *   T2 delta>0, eps=0 (Johnson) — a delta-perturbed exact inclusion into an exact
 *      C* A=M_n -> vtilde an EXACT inclusion (defects machine-zero), c_0*eps=0, and
 *      ||vtilde - v|| <= O(delta).
 *      MUTATION: max_steps=0 -> output not certifiably O(eps) -> abort RED.
 *   T3 eps>0 — a delta-inclusion into a GENUINE eta>0 A (make_mixconj) -> vtilde a
 *      c_0*eps-inclusion; report c_0; every defect <= c_0*eps (by construction of
 *      the certification guard). STAR exercised (eta>0).
 *   T4 THE c_0 UNIVERSALITY CANARY (.tex:484, the headline) — sweep A over
 *      increasing dim (B=M_2..M_5 single block, dim_B 4..25 a 6.25x range; and a
 *      M_2(+)M_2 / M_2(+)M_2(+)M_2 block sweep) and ASSERT the measured c_0 does NOT
 *      grow with dim A/dim B (ratio c0_max/c0_min < 1.6 over the 6.25x range — tight
 *      enough to catch even SUBLINEAR/sqrt-dim growth). A c_0 proportional to dim is
 *      the .tex:484 failure mode (escalate aic-1bc).
 *      MUTATION: inject a dim-factor into c_0 (linear *dim_B OR *sqrt(dim_B)) ->
 *      the ratio assert fires RED (verified, see test_canary comment).
 *   T5 BIJECTIVITY (.tex:1318) — a delta-ISOMORPHISM (dim_B == dim_A) -> vtilde a
 *      c_0*eps-ISOMORPHISM: aic_errreduce_is_bijective returns true and reports
 *      a > 0 (injective + dim match); a NON-square v (dim_B < dim_A) is NOT
 *      bijective.
 *   T6 COLLAPSE-DETECTOR GUARD (F1 blocker) — a delta-hom bounded below on every
 *      basis element but COLLAPSING a combination (v(E_0)=diag(1,0), v(E_1)=rank-1
 *      projection at angle 0.1) passes the OLD basis sweep yet has true inclusion
 *      inf ~0.0998; the SOUND coordinate-matrix sigma_min catches it -> ABORT.
 *      MUTATION: revert the input guard to the basis sweep -> the collapse is NOT
 *      caught (no abort with the delta-inclusion message) -> RED (verified).
 *
 * FULLY CONSTRUCTIVE (no eig): products + certified op-norms, all via dhom.
 */
/* mkstemp / fork / waitpid (T6 stderr-capture fork helper). Must precede all
 * system includes. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include <complex.h>
#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_assoc.h"
#include "aic_dhom.h"
#include "aic_ecstar.h"
#include "aic_errreduce.h"
#include "aic_mat.h"
#include "aic_test.h"
#include "aic_ucp.h"
#include "test_idemp.h"

#undef PREC
#define PREC 80
static double dd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }

__attribute__((unused))
static void suppress_unused_idemp_builders(void)
{
    (void) make_block_cond_exp;
    (void) make_trace_replace;
    (void) make_noiseless_subsystem;
    (void) make_dephasing;
    (void) make_compress_idemp;
}

/* opnorm on ball midpoints (mirrors test_dhom; the certified mid+radius UB lives
 * in the production code). */
static double opnorm_d(const acb_mat_t A, slong prec)
{
    slong r = acb_mat_nrows(A), c = acb_mat_ncols(A);
    acb_mat_t M;
    acb_mat_init(M, r, c);
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++)
            acb_get_mid(acb_mat_entry(M, i, j), acb_mat_entry(A, i, j));
    arb_t e;
    arb_init(e);
    aic_mat_opnorm(e, M, prec);
    double v = dd(e);
    arb_clear(e);
    acb_mat_clear(M);
    return v;
}

/* deterministic pseudo-random complex matrix (no RNG; project rule). */
static void fill_random(acb_mat_t M, unsigned seed)
{
    slong r = acb_mat_nrows(M), c = acb_mat_ncols(M);
    unsigned s = seed * 2654435761u + 1u;
    for (slong i = 0; i < r; i++)
        for (slong j = 0; j < c; j++) {
            s = s * 1103515245u + 12345u;
            double re = ((double) ((s >> 9) & 0xffff) / 32768.0) - 1.0;
            s = s * 1103515245u + 12345u;
            double im = ((double) ((s >> 9) & 0xffff) / 32768.0) - 1.0;
            acb_set_d_d(acb_mat_entry(M, i, j), re, im);
        }
}

/* v(E_i) = U E_i U^dag for B = M_n single block: an exact *-iso M_n -> A=M_n. */
static void set_v_conjugation(aic_dhom_v *v, const acb_mat_t U, slong prec)
{
    slong n = v->n;
    acb_mat_t Ud, Ei, tmp;
    acb_mat_init(Ud, n, n);
    acb_mat_init(Ei, n, n);
    acb_mat_init(tmp, n, n);
    acb_mat_conjugate_transpose(Ud, U);
    for (slong i = 0; i < v->B->dim_B; i++) {
        aic_dhom_B_matunit(Ei, v->B, i);
        acb_mat_mul(tmp, U, Ei, prec);
        acb_mat_mul(v->vE[i], tmp, Ud, prec);
    }
    acb_mat_clear(tmp);
    acb_mat_clear(Ei);
    acb_mat_clear(Ud);
}

/* Project R onto A = span{B_k}: Rp = sum_k <B_k,R>_F B_k (keeps a perturbation IN A). */
static void project_onto_A(acb_mat_t Rp, const aic_ecstar *A, const acb_mat_t R,
                           slong prec)
{
    slong n = A->n;
    acb_mat_zero(Rp);
    acb_t c, z;
    acb_init(c);
    acb_init(z);
    acb_mat_t tmp;
    acb_mat_init(tmp, n, n);
    for (slong k = 0; k < A->dim_A; k++) {
        acb_zero(c);
        for (slong x = 0; x < n; x++)
            for (slong y = 0; y < n; y++) {
                acb_conj(z, acb_mat_entry(A->B[k], x, y));
                acb_mul(z, z, acb_mat_entry(R, x, y), prec);
                acb_add(c, c, z, prec);
            }
        acb_mat_scalar_mul_acb(tmp, A->B[k], c, prec);
        acb_mat_add(Rp, Rp, tmp, prec);
    }
    acb_mat_clear(tmp);
    acb_clear(z);
    acb_clear(c);
}

/* Multi-block conditional expectation (eta=0 idempotent), generalizing
 * make_block_cond_exp to any number of blocks. (Copied from test_dhom.) */
static void make_multiblock_cond_exp(aic_ucp_kraus *phi, const slong *dims,
                                     slong nblk)
{
    slong n = 0;
    for (slong l = 0; l < nblk; l++) n += dims[l];
    aic_ucp_kraus_init(phi, n, n, nblk);
    slong off = 0;
    for (slong l = 0; l < nblk; l++) {
        for (slong i = 0; i < dims[l]; i++)
            acb_set_si(acb_mat_entry(phi->K[l], off + i, off + i), 1);
        off += dims[l];
    }
}

/* Convex mix of make_multiblock_cond_exp with its unitary conjugate: a genuine
 * eta>0 almost-idempotent channel whose A ~ (+)_l M_{dims[l]}. (From test_dhom.) */
static void make_mixconj_bce(aic_ucp_kraus *out, const slong *dims, slong nblk,
                             double t, slong prec)
{
    aic_ucp_kraus base, conj;
    make_multiblock_cond_exp(&base, dims, nblk);
    make_conjugated(&conj, &base, prec);
    slong n = base.dim_H;
    aic_ucp_kraus_init(out, n, n, base.r + conj.r);
    arb_t s;
    arb_init(s);
    arb_set_d(s, sqrt(1.0 - t));
    for (slong a = 0; a < base.r; a++)
        acb_mat_scalar_mul_arb(out->K[a], base.K[a], s, prec);
    arb_set_d(s, sqrt(t));
    for (slong b = 0; b < conj.r; b++)
        acb_mat_scalar_mul_arb(out->K[base.r + b], conj.K[b], s, prec);
    arb_clear(s);
    aic_ucp_kraus_clear(&conj);
    aic_ucp_kraus_clear(&base);
}

/* Build a delta-INCLUSION v: B = (+)_l M_{dims[l]} -> A = Img Phi_tilde of the
 * built eta>0 algebra `ae`. v(E^{(l)}_ab) = Phi_tilde(iota_l(E_ab)) (the natural
 * delta-hom; roff[l] = block l row-offset in M_n), then a *-symmetric delta-
 * perturbation projected into A (mirrors test_dhom run_reduce, so the input is a
 * bona fide delta-inclusion the cor_improvement guard accepts). */
static void build_inclusion(aic_dhom_v *v, aic_dhom_B *B, aic_assoc_ecstar *ae,
                            const slong *dims, slong nblk, const slong *roff,
                            double delta, slong prec)
{
    slong n = ae->A.n;
    aic_dhom_B_init(B, dims, nblk);
    AIC_CHECK_MSG(ae->A.dim_A == B->dim_B, "build_inclusion: dim_A=%ld != dim_B=%ld",
                  (long) ae->A.dim_A, (long) B->dim_B);
    aic_dhom_v_init(v, B, &ae->A);

    acb_mat_t iE, vE;
    acb_mat_init(iE, n, n);
    acb_mat_init(vE, n, n);
    for (slong l = 0; l < nblk; l++) {
        slong dl = dims[l], r0 = roff[l];
        for (slong a = 0; a < dl; a++)
            for (slong b = 0; b < dl; b++) {
                slong i = aic_dhom_B_index(B, l, a, b);
                acb_mat_zero(iE);
                acb_set_si(acb_mat_entry(iE, r0 + a, r0 + b), 1);
                aic_assoc_superop_apply(vE, ae->S_tilde, iE, prec);
                acb_mat_set(v->vE[i], vE);
            }
    }
    /* *-symmetric delta-perturbation, projected into A. */
    arb_t dl_arb, inv, ah;
    arb_init(dl_arb);
    arb_init(inv);
    arb_init(ah);
    arb_set_d(dl_arb, delta);
    arb_set_d(ah, 0.5);
    acb_mat_t R, Rd, Rp;
    acb_mat_init(R, n, n);
    acb_mat_init(Rd, n, n);
    acb_mat_init(Rp, n, n);
    for (slong l = 0; l < nblk; l++)
        for (slong a = 0; a < dims[l]; a++)
            for (slong b = a; b < dims[l]; b++) {
                slong i = aic_dhom_B_index(B, l, a, b);
                slong j = aic_dhom_B_index(B, l, b, a);
                fill_random(R, (unsigned) (201 + i));
                if (a == b) {
                    acb_mat_conjugate_transpose(Rd, R);
                    acb_mat_add(R, R, Rd, prec);
                    acb_mat_scalar_mul_arb(R, R, ah, prec);
                }
                project_onto_A(Rp, &ae->A, R, prec);
                double rn = opnorm_d(Rp, prec);
                arb_set_d(inv, 1.0 / rn);
                acb_mat_scalar_mul_arb(Rp, Rp, inv, prec);
                acb_mat_scalar_mul_arb(Rp, Rp, dl_arb, prec);
                acb_mat_add(v->vE[i], v->vE[i], Rp, prec);
                if (j != i) {
                    acb_mat_conjugate_transpose(Rd, Rp);
                    acb_mat_add(v->vE[j], v->vE[j], Rd, prec);
                }
            }
    acb_mat_clear(Rp);
    acb_mat_clear(Rd);
    acb_mat_clear(R);
    arb_clear(ah);
    arb_clear(inv);
    arb_clear(dl_arb);
    acb_mat_clear(vE);
    acb_mat_clear(iE);
}

/* Single-block convenience: A from make_mixconj(d,m,t), B = M_m at roff 0,
 * run errreduce, return measured c_0 (and report eps / max-defect). */
static double sb_errreduce(slong d, slong m, double t, double delta,
                           double *eps_out, double *maxdef_out, int *bij_out)
{
    aic_ucp_kraus phi;
    make_mixconj(&phi, d, m, t, PREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong dims[1] = {m}, roff[1] = {0};

    aic_dhom_B B;
    aic_dhom_v v;
    build_inclusion(&v, &B, &ae, dims, 1, roff, delta, PREC);

    arb_t epsb;
    arb_init(epsb);
    aic_ecstar_defect_assoc(epsb, &ae.A, PREC);
    double eps = dd(epsb);

    aic_dhom_v vt;
    aic_dhom_v_init(&vt, &B, &ae.A);
    aic_errreduce_defects defs;
    aic_errreduce_defects_init(&defs);
    arb_t c0;
    arb_init(c0);
    aic_errreduce(&vt, c0, &defs, &v, eps, /*tol_abs*/ 1e-13, /*max_steps*/ 30, PREC);

    arb_t maxd;
    arb_init(maxd);
    aic_errreduce_defects_max(maxd, &defs);
    if (eps_out) *eps_out = eps;
    if (maxdef_out) *maxdef_out = dd(maxd);
    if (bij_out) *bij_out = aic_errreduce_is_bijective(NULL, &vt, PREC);
    double c0v = dd(c0);

    arb_clear(maxd);
    arb_clear(c0);
    aic_errreduce_defects_clear(&defs);
    aic_dhom_v_clear(&vt);
    arb_clear(epsb);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return c0v;
}

/* Multi-block convenience: A from make_mixconj_bce(dims), B = (+)_l M_{dims[l]}. */
static double mb_errreduce(const slong *dims, slong nblk, double t, double delta,
                           double *eps_out)
{
    aic_ucp_kraus phi;
    make_mixconj_bce(&phi, dims, nblk, t, PREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong roff[8];
    AIC_CHECK_MSG(nblk <= 8, "mb_errreduce: too many blocks");
    slong o = 0;
    for (slong l = 0; l < nblk; l++) { roff[l] = o; o += dims[l]; }

    aic_dhom_B B;
    aic_dhom_v v;
    build_inclusion(&v, &B, &ae, dims, nblk, roff, delta, PREC);

    arb_t epsb;
    arb_init(epsb);
    aic_ecstar_defect_assoc(epsb, &ae.A, PREC);
    double eps = dd(epsb);

    aic_dhom_v vt;
    aic_dhom_v_init(&vt, &B, &ae.A);
    arb_t c0;
    arb_init(c0);
    aic_errreduce(&vt, c0, NULL, &v, eps, 1e-13, 30, PREC);
    double c0v = dd(c0);
    if (eps_out) *eps_out = eps;

    arb_clear(c0);
    aic_dhom_v_clear(&vt);
    arb_clear(epsb);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
    return c0v;
}

/* Run `child` in a forked process with stderr captured to a temp file; after it
 * exits, assert it died by SIGABRT and the captured stderr CONTAINS `expect`
 * (mirrors test_projection.c expect_abort_with_msg; the message-grep gives the
 * fail-loud test teeth — a DIFFERENT guard firing changes the message -> RED). */
static void expect_abort_with_msg(void (*child)(void), const char *expect,
                                  const char *name)
{
    char tmpl[] = "/tmp/aic_errreduce_failXXXXXX";
    int fd = mkstemp(tmpl);
    AIC_CHECK_MSG(fd >= 0, "%s: mkstemp failed", name);
    pid_t pid = fork();
    AIC_CHECK_MSG(pid >= 0, "%s: fork failed", name);
    if (pid == 0) {
        dup2(fd, STDERR_FILENO);
        close(fd);
        child();
        _exit(0);   /* reached only if NO abort fired (a test failure) */
    }
    int status = 0;
    waitpid(pid, &status, 0);
    close(fd);
    int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
    char buf[4096];
    size_t got = 0;
    FILE *f = fopen(tmpl, "r");
    if (f) { got = fread(buf, 1, sizeof(buf) - 1, f); fclose(f); }
    buf[got] = '\0';
    unlink(tmpl);
    int has_msg = (strstr(buf, expect) != NULL);
    printf("%s: child %s; stderr %s expected substring \"%s\"\n", name,
           aborted ? "SIGABRT" : "did NOT abort (WRONG)",
           has_msg ? "CONTAINS" : "MISSING", expect);
    AIC_CHECK_MSG(aborted, "%s: child did not abort via SIGABRT", name);
    AIC_CHECK_MSG(has_msg, "%s: captured stderr does NOT contain \"%s\" (got: "
                  "%.200s) -- the INTENDED guard did not fire", name, expect, buf);
}

/* Build the reviewer's COLLAPSING-v fixture (the F1 BLOCKER witness).
 *   B = C (+) C  (dims {1,1}, dim_B = 2),  A = M_2 (identity channel, dim_A = 4).
 *   v(E_0) = diag(1,0),   v(E_1) = |u><u| with u = (cos t, sin t), t = 0.1.
 * BASIS sweep: ||v(E_0)||_op = ||v(E_1)||_op = 1 (rank-1 projection), so the basis
 * lower bound a = 1.0 PASSES the old basis-sweep guard. But v COLLAPSES the
 * direction E_0 - E_1: ||v(E_0 - E_1)||_op = ||diag(1,0) - |u><u|||_op = |sin t| =
 * 0.0998, while ||E_0 - E_1||_op = 1, so the TRUE inclusion infimum <= 0.0998 <<
 * 0.5 — the delta-inclusion hypothesis (1-delta)||X|| <= ||v(X)|| is VIOLATED. The
 * coordinate-matrix sigma_min SEES this (sigma_min(M) ~ 0.0998), so the SOUND guard
 * aborts; the old basis-sweep guard did NOT (the BLOCKER). Both v(E_i) are in
 * A = M_2 (full algebra), so this is a legitimate delta-hom whose ONLY defect is
 * the collapse. */
static void build_collapsing_v(aic_dhom_v *v, aic_dhom_B *B, aic_assoc_ecstar *ae)
{
    slong dims[2] = {1, 1};
    aic_dhom_B_init(B, dims, 2);          /* B = C (+) C, dim_B = 2 */
    aic_dhom_v_init(v, B, &ae->A);
    double t = 0.1, c = cos(t), s = sin(t);
    /* v(E_0) = diag(1, 0). */
    acb_mat_zero(v->vE[0]);
    acb_set_si(acb_mat_entry(v->vE[0], 0, 0), 1);
    /* v(E_1) = |u><u|, u = (c, s):  [[c^2, c s], [c s, s^2]]. */
    acb_mat_zero(v->vE[1]);
    acb_set_d(acb_mat_entry(v->vE[1], 0, 0), c * c);
    acb_set_d(acb_mat_entry(v->vE[1], 0, 1), c * s);
    acb_set_d(acb_mat_entry(v->vE[1], 1, 0), c * s);
    acb_set_d(acb_mat_entry(v->vE[1], 1, 1), s * s);
}

static void t6_collapse_child(void)
{
    aic_ucp_kraus phi;
    make_identity(&phi, 2);               /* A = M_2 (eta=0, full algebra) */
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    aic_dhom_B B;
    aic_dhom_v v;
    build_collapsing_v(&v, &B, &ae);
    aic_dhom_v vt;
    aic_dhom_v_init(&vt, &B, &ae.A);
    arb_t c0;
    arb_init(c0);
    /* MUST abort at the input delta-inclusion guard (sigma_min ~ 0.0998 < 0.5). */
    aic_errreduce(&vt, c0, NULL, &v, /*eps*/ 0.0, /*tol_abs*/ 1e-13,
                  /*max_steps*/ 8, PREC);
    _exit(0);   /* reached only if the (unsound) guard let the collapse through */
}

/* ====================================================================== T6 */
/* THE F1 BLOCKER FIX: the input delta-inclusion guard must be SOUND (sigma_min of
 * the coordinate matrix), not the basis-blind basis sweep. A v that is bounded
 * below on every basis element but COLLAPSES a general combination (v(E_0 - E_1))
 * must ABORT. Also: a GENUINE delta-inclusion (T1/T3/T5 fixtures) must still pass
 * (no false abort) — re-checked there. */
static void test_collapse_guard(void)
{
    printf("-- T6 collapse-detector guard (F1 blocker) --\n");
    /* Diagnostic: report the basis-sweep a (passes) vs the sigma_min (catches). */
    aic_ucp_kraus phi;
    make_identity(&phi, 2);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    aic_dhom_B B;
    aic_dhom_v v;
    build_collapsing_v(&v, &B, &ae);
    arb_t a_basis, smin;
    arb_init(a_basis);
    arb_init(smin);
    aic_dhom_prop_bounds(NULL, a_basis, NULL, &v, PREC);
    aic_dhom_v_sigma_min(smin, &v, PREC);
    printf("   collapsing v: basis-sweep a = %.6f (PASSES 0.5 guard), "
           "sigma_min = %.6f (CATCHES)\n", dd(a_basis), dd(smin));
    AIC_CHECK_MSG(dd(a_basis) >= 0.5, "T6: basis-sweep a=%.6f should be >= 0.5 "
                  "(the witness must defeat the OLD guard)", dd(a_basis));
    AIC_CHECK_MSG(dd(smin) < 0.5, "T6: sigma_min=%.6f should be < 0.5 (the SOUND "
                  "guard must catch the collapse)", dd(smin));
    AIC_CHECK_MSG(fabs(dd(smin) - sin(0.1)) < 1e-6, "T6: sigma_min=%.6f != "
                  "|sin 0.1|=%.6f (expected true inclusion inf)", dd(smin), sin(0.1));
    arb_clear(smin);
    arb_clear(a_basis);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);

    /* the guard must ABORT on this fixture (fail-loud, message grep). */
    expect_abort_with_msg(t6_collapse_child, "NOT a delta-inclusion",
                          "T6 collapse-abort");
}

/* ====================================================================== T1 */
/* eta=0 ORACLE: A = M_n identity channel, v exact *-iso (a 0-inclusion).
 * errreduce returns vtilde = v unchanged, all defects machine-zero, c_0*eps = 0. */
static void test_eta0_oracle(void)
{
    printf("-- T1 eta=0 oracle --\n");
    for (slong n = 2; n <= 3; n++) {
        aic_ucp_kraus phi;
        make_identity(&phi, n);
        aic_assoc_ecstar ae;
        aic_assoc_ecstar_from_phi(&ae, &phi, PREC);

        slong dims[1] = {n};
        aic_dhom_B B;
        aic_dhom_B_init(&B, dims, 1);
        aic_dhom_v v;
        aic_dhom_v_init(&v, &B, &ae.A);
        acb_mat_t U;
        acb_mat_init(U, n, n);
        make_fixed_unitary(U, n, PREC);
        set_v_conjugation(&v, U, PREC);

        aic_dhom_v vt;
        aic_dhom_v_init(&vt, &B, &ae.A);
        aic_errreduce_defects defs;
        aic_errreduce_defects_init(&defs);
        arb_t c0;
        arb_init(c0);
        aic_errreduce(&vt, c0, &defs, &v, /*eps*/ 0.0, /*tol_abs*/ 1e-15,
                      /*max_steps*/ 8, PREC);

        AIC_CHECK_MSG(dd(c0) < 1e-14, "T1(n=%ld): c_0 %.3e != 0 (eta=0)",
                      (long) n, dd(c0));
        arb_t maxd;
        arb_init(maxd);
        aic_errreduce_defects_max(maxd, &defs);
        AIC_CHECK_MSG(dd(maxd) < 1e-14, "T1(n=%ld): max inclusion-defect %.3e != 0",
                      (long) n, dd(maxd));
        printf("   n=%ld: c_0=%.3e  max-defect=%.3e (mult=%.2e nrm=%.2e low=%.2e "
               "unit=%.2e)\n", (long) n, dd(c0), dd(maxd), dd(defs.mult_def),
               dd(defs.norm_excess), dd(defs.lower_gap), dd(defs.unit_def));

        /* vtilde unchanged from v (exact hom: lem_approx is a no-op). */
        double maxchg = 0.0;
        acb_mat_t diff;
        acb_mat_init(diff, n, n);
        for (slong i = 0; i < B.dim_B; i++) {
            acb_mat_sub(diff, vt.vE[i], v.vE[i], PREC);
            double c = opnorm_d(diff, PREC);
            if (c > maxchg) maxchg = c;
        }
        AIC_CHECK_MSG(maxchg < 1e-15, "T1(n=%ld): vtilde changed v by %.3e",
                      (long) n, maxchg);
        acb_mat_clear(diff);

        /* MUTATION (verified RED 2026-05-30): collapse a basis op of v (vE[1]=0),
         * making v a NON-inclusion (basis-sweep lower bound a -> 0). Then
         * aic_errreduce aborts at the input-delta-inclusion guard. Re-enable to
         * confirm:
         *   acb_mat_zero(v.vE[1]);
         *   aic_errreduce(&vt, c0, &defs, &v, 0.0, 1e-15, 8, PREC);  // SIGABRT
         * Restored: the guard fires "input is NOT a delta-inclusion (a < 0.5)". */

        arb_clear(maxd);
        arb_clear(c0);
        aic_errreduce_defects_clear(&defs);
        acb_mat_clear(U);
        aic_dhom_v_clear(&vt);
        aic_dhom_v_clear(&v);
        aic_dhom_B_clear(&B);
        aic_assoc_ecstar_clear(&ae);
        aic_ucp_kraus_clear(&phi);
    }
}

/* ====================================================================== T2 */
/* delta>0, eps=0 (Johnson): a delta-perturbed exact inclusion into an exact C*
 * A=M_n -> vtilde an EXACT inclusion (defects machine-zero), c_0*eps=0, and
 * ||vtilde - v_pert|| <= O(delta). */
static void test_delta_eps0(void)
{
    printf("-- T2 delta>0, eps=0 (Johnson) --\n");
    slong n = 2;
    double delta = 0.02;
    aic_ucp_kraus phi;
    make_identity(&phi, n);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong dims[1] = {n};
    aic_dhom_B B;
    aic_dhom_B_init(&B, dims, 1);
    aic_dhom_v v;
    aic_dhom_v_init(&v, &B, &ae.A);
    acb_mat_t U;
    acb_mat_init(U, n, n);
    make_fixed_unitary(U, n, PREC);
    set_v_conjugation(&v, U, PREC);

    /* *-symmetric delta-perturbation (keeps v *-linear), into the exact A=M_n. */
    arb_t dlt, inv, ah;
    arb_init(dlt);
    arb_init(inv);
    arb_init(ah);
    arb_set_d(dlt, delta);
    arb_set_d(ah, 0.5);
    acb_mat_t R, Rn, Rd;
    acb_mat_init(R, n, n);
    acb_mat_init(Rn, n, n);
    acb_mat_init(Rd, n, n);
    for (slong a = 0; a < n; a++)
        for (slong b = 0; b < n; b++) {
            slong i = a * n + b, j = b * n + a;
            if (j < i) continue;
            fill_random(R, (unsigned) (101 + i));
            if (a == b) {
                acb_mat_conjugate_transpose(Rd, R);
                acb_mat_add(R, R, Rd, PREC);
                acb_mat_scalar_mul_arb(R, R, ah, PREC);
            }
            double rn = opnorm_d(R, PREC);
            arb_set_d(inv, 1.0 / rn);
            acb_mat_scalar_mul_arb(Rn, R, inv, PREC);
            acb_mat_scalar_mul_arb(Rn, Rn, dlt, PREC);
            acb_mat_add(v.vE[i], v.vE[i], Rn, PREC);
            if (j != i) {
                acb_mat_conjugate_transpose(Rd, Rn);
                acb_mat_add(v.vE[j], v.vE[j], Rd, PREC);
            }
        }

    aic_dhom_v vt;
    aic_dhom_v_init(&vt, &B, &ae.A);
    aic_errreduce_defects defs;
    aic_errreduce_defects_init(&defs);
    arb_t c0;
    arb_init(c0);
    /* eps=0 (exact C* A): errreduce drives the defects to machine-zero. */
    aic_errreduce(&vt, c0, &defs, &v, /*eps*/ 0.0, /*tol_abs*/ 1e-13,
                  /*max_steps*/ 40, PREC);
    arb_t maxd;
    arb_init(maxd);
    aic_errreduce_defects_max(maxd, &defs);
    printf("   delta=%.2f: c_0=%.3e  max-defect=%.3e\n", delta, dd(c0), dd(maxd));
    AIC_CHECK_MSG(dd(c0) < 1e-12, "T2: c_0 %.3e != 0 (eps=0)", dd(c0));
    AIC_CHECK_MSG(dd(maxd) < 1e-12, "T2: exact-inclusion max-defect %.3e != 0 "
                  "(Johnson eps=0 limit)", dd(maxd));

    /* ||vtilde - v_pert|| <= O(delta) */
    double maxchg = 0.0;
    acb_mat_t diff;
    acb_mat_init(diff, n, n);
    for (slong i = 0; i < B.dim_B; i++) {
        acb_mat_sub(diff, vt.vE[i], v.vE[i], PREC);
        double c = opnorm_d(diff, PREC);
        if (c > maxchg) maxchg = c;
    }
    printf("   ||vtilde - v_pert||_max = %.3e = %.2f * delta\n", maxchg,
           maxchg / delta);
    AIC_CHECK_MSG(maxchg <= 10.0 * delta, "T2: ||vtilde - v|| %.3e not O(delta)",
                  maxchg);
    acb_mat_clear(diff);

    /* MUTATION (verified RED): max_steps=0 leaves vtilde=v_pert whose mult-defect
     * = delta >> 0; the certification guard (max-defect <= C0_CERT*tol) aborts.
     * Re-enable to confirm:
     *   aic_errreduce(&vt2, c0, &defs, &v, 0.0, 1e-13, 0, PREC);  // SIGABRT
     * Restored. */

    arb_clear(maxd);
    arb_clear(c0);
    aic_errreduce_defects_clear(&defs);
    aic_dhom_v_clear(&vt);
    acb_mat_clear(Rd);
    acb_mat_clear(Rn);
    acb_mat_clear(R);
    arb_clear(ah);
    arb_clear(inv);
    arb_clear(dlt);
    acb_mat_clear(U);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

/* ====================================================================== T3 */
/* eps>0 GENUINE almost-idempotent A: a delta-inclusion into make_mixconj(5,3) ->
 * vtilde a c_0*eps-inclusion; report c_0; defects <= c_0*eps (the certification
 * guard enforces it). STAR exercised (eta>0). */
static void test_eps_pos(void)
{
    printf("-- T3 eps>0 genuine almost-idempotent A --\n");
    double eps, maxdef;
    int bij;
    double c0 = sb_errreduce(/*d*/ 5, /*m*/ 3, /*t*/ 0.015, /*delta*/ 0.05, &eps,
                             &maxdef, &bij);
    printf("   eps_assoc=%.3e  c_0=%.3f  max-defect=%.3e (= %.2f * eps)\n", eps,
           c0, maxdef, maxdef / eps);
    AIC_CHECK_MSG(c0 > 0.0, "T3: c_0 = %.3f not positive (eps>0)", c0);
    AIC_CHECK_MSG(maxdef <= AIC_ERRREDUCE_C0_CERT * eps, "T3: max-defect %.3e > "
                  "C0_CERT*eps (not an O(eps)-inclusion)", maxdef);
    /* the defect = c_0 * eps by construction; confirm the round-trip. */
    AIC_CHECK_MSG(fabs(maxdef - c0 * eps) < 1e-12, "T3: maxdef %.3e != c_0*eps "
                  "%.3e", maxdef, c0 * eps);
}

/* ====================================================================== T4 */
/* THE c_0 UNIVERSALITY CANARY (.tex:484/461). The measured c_0 must NOT grow with
 * dim A / dim B. Two sweeps: (A) single-block B=M_2..M_5 (dim_B 4..25, a 6.25x
 * range); (B) number-of-blocks m=1,2,3 of M_2. Assert the c_0 ratio across each
 * sweep is bounded by a constant. A c_0 proportional to dim is the .tex:484 failure
 * mode the paper warns against (escalate aic-1bc).
 *
 * THE RATIO THRESHOLD HAS TEETH AGAINST SUBLINEAR GROWTH (F2 fix). The observed
 * (A) spread is small (raw c_0 DECREASES with dim: 2.71 -> 2.22 -> 2.07 -> 1.96
 * over the sweep, ratio ~1.38), so the threshold is tightened to 1.6 over the
 * WIDENED 6.25x dim range. This catches even sqrt(dim) growth: an injected
 * c_0 *= sqrt(dim_B) inflates the M_5/M_2 ratio by sqrt(25/4) = 2.5x on top of the
 * raw spread, far above 1.6 (the OLD 1.525-over-4x-range escape no longer hides
 * sublinear growth).
 *
 * MUTATION (verified RED 2026-05-30): inject a dim-factor into the measured c_0 in
 * aic_errreduce, after computing c0_out:
 *     arb_mul_si(c0_out, c0_out, vtilde_out->B->dim_B, prec);          (linear)  or
 *     { arb_t s; arb_init(s); arb_sqrt_ui(s, (ulong)vtilde_out->B->dim_B, prec);
 *       arb_mul(c0_out, c0_out, s, prec); arb_clear(s); }              (sqrt-dim)
 * The sqrt-dim injection drives the (A) ratio to 1.807 (= sqrt(25/4) * 1.384/1.92,
 * measured) > 1.6 -> the (A) assert fires RED (T4(A) line, INDEPENDENTLY of T3 — the
 * T3 round-trip maxdef==c0*eps also trips on a c_0 tamper, but disabling T3 confirms
 * T4(A) is RED-able on its own). Linear injection RED a fortiori. Restored. This
 * proves the canary catches a dimension-dependent constant, linear OR sublinear. */
static void test_canary(void)
{
    printf("-- T4 c_0 universality canary (.tex:484) --\n");
    double eps;

    /* (A) single-block dimension sweep B = M_2..M_5 (dim_B 4,9,16,25 — a 6.25x
     * range; sb_errreduce(d=m+2, m) — a fixed compression overhead d-m=2 across the
     * sweep, so each A = M_m, dim_A = m^2 = dim_B and the ONLY varying quantity is
     * the dim, isolating the .tex:484 dimension dependence. M_6 (n=8) is dropped to
     * keep the PREC=80 run in the seconds regime; 6.25x is still >= the 6x the F2
     * review asked for, enough teeth against sqrt-dim growth). */
    const int NSWEEP = 4;
    double Cdim[4];
    slong msweep[4] = {2, 3, 4, 5};
    for (int i = 0; i < NSWEEP; i++) {
        slong m = msweep[i];
        Cdim[i] = sb_errreduce(m + 2, m, 0.015, 0.05, &eps, NULL, NULL);
        printf("   (A) B=M_%ld  dim_B=%ld: eps=%.3e c_0=%.3f\n", (long) m,
               (long) (m * m), eps, Cdim[i]);
    }
    double dmin = Cdim[0], dmax = Cdim[0];
    for (int i = 1; i < NSWEEP; i++) {
        if (Cdim[i] < dmin) dmin = Cdim[i];
        if (Cdim[i] > dmax) dmax = Cdim[i];
    }
    printf("   (A) c_0 in [%.3f, %.3f], ratio c0_max/c0_min = %.3f (dim_B 4..25)\n",
           dmin, dmax, dmax / dmin);
    AIC_CHECK_MSG(dmax / dmin < 1.6, "T4(A): c_0 ratio %.3f grows with dim_B (the "
                  ".tex:484 failure mode — even SUBLINEAR/sqrt-dim growth; the c_0 "
                  "must be dimension-independent — escalate aic-1bc)", dmax / dmin);

    /* (B) number-of-blocks sweep m = 1,2,3 of M_2. The embedded-sum diagonal has
     * ||D||_proj = m (FINDINGS A2/D2), so a worst-case w' bound is O(m*delta); we
     * assert the achieved c_0 does NOT carry that m-factor as growth. */
    double Cm[3];
    Cm[0] = sb_errreduce(4, 2, 0.015, 0.05, &eps, NULL, NULL); /* m=1 */
    printf("   (B) m=1 B=M_2          dim_B=4 : eps=%.3e c_0=%.3f\n", eps, Cm[0]);
    {
        slong d2[2] = {2, 2};
        Cm[1] = mb_errreduce(d2, 2, 0.015, 0.05, &eps);
        printf("   (B) m=2 B=M_2(+)M_2    dim_B=8 : eps=%.3e c_0=%.3f\n", eps, Cm[1]);
        slong d3[3] = {2, 2, 2};
        Cm[2] = mb_errreduce(d3, 3, 0.015, 0.05, &eps);
        printf("   (B) m=3 B=M_2(+)M_2(+)M_2 dim_B=12: eps=%.3e c_0=%.3f\n", eps,
               Cm[2]);
    }
    double cmmax = Cm[0];
    for (int i = 1; i < 3; i++) if (Cm[i] > cmmax) cmmax = Cm[i];
    printf("   (B) max_m c_0 = %.3f; c_0(1)=%.3f c_0(2)=%.3f c_0(3)=%.3f -> c_0 %s "
           "with m\n", cmmax, Cm[0], Cm[1], Cm[2],
           Cm[2] > Cm[0] ? "GROWS" : "does NOT grow");
    AIC_CHECK_MSG(cmmax < 6.0, "T4(B): max_m c_0 = %.3f exceeds block-count-"
                  "independent ceiling (the O(m*delta) w' bound manifesting)", cmmax);
    AIC_CHECK_MSG(Cm[2] <= 1.3 * Cm[0], "T4(B): c_0(m=3)=%.3f > c_0(m=1)=%.3f — "
                  "c_0 GROWS with the number of blocks (escalate aic-1bc)", Cm[2],
                  Cm[0]);
}

/* ====================================================================== T5 */
/* BIJECTIVITY (.tex:1318). A delta-ISOMORPHISM (dim_B == dim_A) -> vtilde a
 * c_0*eps-ISOMORPHISM: aic_errreduce_is_bijective returns true and a > 0
 * (injective + dim match). A non-square inclusion (dim_B < dim_A) is NOT bijective. */
static void test_bijectivity(void)
{
    printf("-- T5 bijectivity --\n");

    /* (a) square (dim_B == dim_A): the make_mixconj(d,m) A has dim_A = m^2 and
     * B = M_m has dim_B = m^2, so the natural inclusion is a delta-ISOMORPHISM.
     * sb_errreduce reports the bijectivity flag of vtilde. */
    double eps, maxdef;
    int bij;
    double c0 = sb_errreduce(5, 3, 0.015, 0.05, &eps, &maxdef, &bij);
    printf("   square M_3 -> A (dim_B=dim_A=9): bijective=%d  c_0=%.3f  eps=%.3e\n",
           bij, c0, eps);
    AIC_CHECK_MSG(bij == 1, "T5: square delta-isomorphism vtilde NOT reported "
                  "bijective (injective + dim match should hold)");

    /* (b) explicitly confirm the injective + dim-match decomposition on a fresh
     * vtilde, exposing the certified lower bound a > 0. */
    aic_ucp_kraus phi;
    make_mixconj(&phi, 5, 3, 0.015, PREC);
    aic_assoc_ecstar ae;
    aic_assoc_ecstar_from_phi(&ae, &phi, PREC);
    slong dims[1] = {3}, roff[1] = {0};
    aic_dhom_B B;
    aic_dhom_v v;
    build_inclusion(&v, &B, &ae, dims, 1, roff, 0.05, PREC);
    arb_t epsb;
    arb_init(epsb);
    aic_ecstar_defect_assoc(epsb, &ae.A, PREC);
    aic_dhom_v vt;
    aic_dhom_v_init(&vt, &B, &ae.A);
    arb_t c0b;
    arb_init(c0b);
    aic_errreduce(&vt, c0b, NULL, &v, dd(epsb), 1e-13, 30, PREC);

    arb_t a;
    arb_init(a);
    int b2 = aic_errreduce_is_bijective(a, &vt, PREC);
    printf("   is_bijective: %d  (lower bound a = %.6f > 0, dim_B=%ld dim_A=%ld)\n",
           b2, dd(a), (long) B.dim_B, (long) ae.A.dim_A);
    AIC_CHECK_MSG(b2 == 1, "T5(b): vtilde not bijective");
    AIC_CHECK_MSG(dd(a) > 0.5, "T5(b): lower bound a=%.6f not > 0.5 (not injective)",
                  dd(a));
    AIC_CHECK_MSG(B.dim_B == ae.A.dim_A, "T5(b): dim mismatch %ld != %ld",
                  (long) B.dim_B, (long) ae.A.dim_A);
    /* the inverse of a delta-iso is a (delta+O(delta^2))-iso (.tex:458); since
     * vtilde is a c_0*eps-iso with a >= 1-O(eps) > 0, the inverse exists and is
     * bounded by 1/a = 1+O(eps). Confirm a is in the O(eps)-isometry band. */
    AIC_CHECK_MSG(dd(a) >= 1.0 - AIC_ERRREDUCE_C0_CERT * dd(epsb),
                  "T5(b): a=%.6f below 1-O(eps) band (inverse not O(eps)-iso)", dd(a));

    /* (c) NON-SQUARE negative case: an INJECTIVE map B' = M_2 (dim_B=4) into the
     * SAME A (dim_A=9). v'(E_ab) = Phi_tilde(iota(E_ab)) for the top-left 2x2 block
     * is bounded below (a > 0), but dim_B != dim_A, so is_bijective must return
     * FALSE (the dim-match veto). This gives the dim_match condition teeth: a
     * MUTATION dropping `&& dim_match` from aic_errreduce_is_bijective would make
     * this case wrongly report bijective -> RED (verified 2026-05-30). */
    slong dims2[1] = {2};
    aic_dhom_B B2;
    aic_dhom_B_init(&B2, dims2, 1);
    aic_dhom_v vns;
    aic_dhom_v_init(&vns, &B2, &ae.A);
    {
        slong nn = ae.A.n;
        acb_mat_t iE, vE;
        acb_mat_init(iE, nn, nn);
        acb_mat_init(vE, nn, nn);
        for (slong aa = 0; aa < 2; aa++)
            for (slong bb = 0; bb < 2; bb++) {
                slong i = aic_dhom_B_index(&B2, 0, aa, bb);
                acb_mat_zero(iE);
                acb_set_si(acb_mat_entry(iE, aa, bb), 1); /* iota(E_ab), top-left */
                aic_assoc_superop_apply(vE, ae.S_tilde, iE, PREC);
                acb_mat_set(vns.vE[i], vE);
            }
        acb_mat_clear(vE);
        acb_mat_clear(iE);
    }
    arb_t ans;
    arb_init(ans);
    int bns = aic_errreduce_is_bijective(ans, &vns, PREC);
    printf("   non-square M_2 -> A (dim_B=4 != dim_A=%ld): bijective=%d (a=%.6f)\n",
           (long) ae.A.dim_A, bns, dd(ans));
    AIC_CHECK_MSG(bns == 0, "T5(c): non-square (dim_B=4 != dim_A=9) wrongly "
                  "reported bijective (the dim-match veto failed)");
    AIC_CHECK_MSG(dd(ans) > 0.5, "T5(c): non-square map a=%.6f not injective "
                  "(should be bounded below; the veto must be the DIM check, not "
                  "injectivity)", dd(ans));
    arb_clear(ans);
    aic_dhom_v_clear(&vns);
    aic_dhom_B_clear(&B2);

    arb_clear(a);
    arb_clear(c0b);
    aic_dhom_v_clear(&vt);
    arb_clear(epsb);
    aic_dhom_v_clear(&v);
    aic_dhom_B_clear(&B);
    aic_assoc_ecstar_clear(&ae);
    aic_ucp_kraus_clear(&phi);
}

int main(void)
{
    flint_set_num_threads(1);
    test_eta0_oracle();
    test_delta_eps0();
    test_eps_pos();
    test_canary();
    test_bijectivity();
    test_collapse_guard();
    aic_test_report("test_errreduce");
    return 0;
}
