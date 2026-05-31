/* aic_cstar_internal.h — INTERNAL shared declarations for the I5 master loop
 * (bead aic-097), split across src/aic_cstar_build.c (driver + Stage 1) and
 * src/aic_cstar_stages.c (Stage 2 + Stage 3) under the 200-LOC/file rule (Rule 10).
 * NOT part of the public API (include/aic_cstar.h is); these are the cross-file
 * helpers the two source files share. See docs/research/cstar_masterloop_spec.md.
 */
#ifndef AIC_CSTAR_INTERNAL_H
#define AIC_CSTAR_INTERNAL_H

#include <stdio.h>
#include <stdlib.h>

#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic_dhom.h"
#include "aic_ecstar.h"

#ifdef __cplusplus
extern "C" {
#endif

/* eps' = O(eps) factor (spec §4.4): a dimension-independent constant giving a
 * margin over the empirical ~2.7*eps wrapper-defect upper bound from errreduce
 * (FINDINGS §D2). The S_P sub-algebras are eps'-C* algebras for any c_0*eps
 * projection P; this is the errreduce FLOOR passed in Stage 2/3 (NOT bare eps). */
#define AIC_STAGE2_EPS_PRIME_FACTOR 4.0

/* Fail-loud helper: abort with a printf-style message if cond holds (Rule 4).
 * Used for the contradiction/canary/nontriviality guards that are stop conditions
 * (not ordinary assert()s — these print diagnostic numbers). */
#define AIC_FAIL_IF(cond, ...)                                                 \
    do {                                                                       \
        if (cond) {                                                            \
            fprintf(stderr, "aic_cstar_build FATAL at %s:%d: ", __FILE__,      \
                    __LINE__);                                                 \
            fprintf(stderr, __VA_ARGS__);                                      \
            fprintf(stderr, "\n");                                             \
            abort();                                                           \
        }                                                                      \
    } while (0)

/* The Stage-1 working set of projections {P_1,...,P_m} (each n x n, in A). A
 * growable array; split() replaces P_m by the pair {P', P''}. */
typedef struct {
    acb_mat_t *P;   /* m operators, each n x n */
    slong m;        /* current count */
    slong cap;      /* allocated capacity */
    slong n;        /* ambient dim (A->n) */
} aic_cstar_plist;

void aic_cstar_plist_init(aic_cstar_plist *list, slong n);
void aic_cstar_plist_clear(aic_cstar_plist *list);
void aic_cstar_plist_push(aic_cstar_plist *list, const acb_mat_t P);

/* The c_0 universality gate (FINDINGS §D2): abort if c0_this > 3*c0_nominal + 5. */
void aic_cstar_c0_gate(double c0_this, double c0_nominal, const char *where);

/* Stage 1 (src/aic_cstar_build.c): fill `list` with the maximal one-dimensional
 * commutative skeleton, write the nominal c_0 (first cor_improvement). */
void aic_cstar_stage1_greedy(aic_cstar_plist *list, double *c0_nominal,
                             const aic_ecstar *A, double eps, double tol_abs,
                             slong max_steps, slong prec);

/* Equivalence-class partition (src/aic_cstar_classes.c, .tex:1428): union-find on
 * the m projections via dim S_{P_j,P_k} == 1. Writes class_id[0..m-1] (a class
 * index in [0,*num_classes)) and *num_classes. class_id must be m-long. ASSERTS
 * (Stage-3 precondition) S_{P_C,P_D}=0 for distinct classes. */
void aic_cstar_compute_classes(slong *class_id, slong *num_classes,
                               const aic_ecstar *A, const aic_cstar_plist *list,
                               slong prec);

/* cor_improvement for a map into an S_P sub-algebra (src/aic_cstar_classes.c):
 * certifies the four inclusion-defects with unit_def against the SUPPLIED unit
 * (Ptilde for a wrapper, the running P_total for a Stage-3 merge), NOT the ambient
 * 1_n (FINDINGS §C7). vt_out init'd over v_in's B and A; c0_out caller-init'd. */
void aic_cstar_errreduce_unit(aic_dhom_v *vt_out, arb_t c0_out,
                              const aic_dhom_v *v_in, const acb_mat_t unit,
                              double eps, double tol_abs, slong max_steps, slong prec);

/* Stages 2+3 (src/aic_cstar_stages.c): per-class inductive extension then merge.
 * Fills B_out, v_out (OWNED, allocated here) and iso_def (NULL to skip). */
void aic_cstar_stages_run(aic_dhom_B *B_out, aic_dhom_v *v_out, arb_t iso_def,
                          const aic_ecstar *A, const aic_cstar_plist *list,
                          const slong *class_id, slong num_classes,
                          double c0_nominal, double eps, double eps_prime,
                          double tol_abs, slong max_steps, slong prec);

#ifdef __cplusplus
}
#endif

#endif /* AIC_CSTAR_INTERNAL_H */
