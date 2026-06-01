/* aic_channels_unitary.c — mixed-unitary public channel constructors (bead
 * aic-7hg): mix_unitaries (the general unital channel), group_twirl (uniform
 * weights), and the qubit pauli channel. See include/aic/aic_channels.h for the
 * Heisenberg convention.
 *
 * THE DAG-ORDER (load-bearing). The desired OBSERVABLE action is
 *     Phi(X) = sum_j p_j U_j X U_j^dag.
 * In the stored form Phi(X) = sum_j K_j^dag X K_j this requires K_j^dag = U_j,
 * i.e. K_j = sqrt(p_j) U_j^dag. Then sum_j K_j^dag K_j = sum_j p_j U_j U_j^dag =
 * sum_j p_j 1 = 1 (unital). Storing K_j = sqrt(p_j) U_j (the WRONG order) would
 * give the TRANSPOSED channel Phi(X)=sum p_j U_j^dag X U_j; the unitarity +
 * sum-to-1 validation does not catch that, so the order is pinned here and
 * cross-checked in tests (a closed-group twirl is idempotent only with the
 * correct convention).
 *
 * pauli is the d=2 special case: sigma_k are Hermitian (U_k^dag = U_k = sigma_k),
 * so K_k = sqrt(p_k) sigma_k and the dag-order is moot; built directly without
 * going through mix_unitaries to avoid allocating the four 2x2 acb_mat unitaries.
 */
#include <math.h>

#include <flint/acb.h>
#include <flint/acb_mat.h>
#include <flint/arb.h>

#include "aic/aic_channels.h"
#include "aic/aic_ucp.h"
#include "aic_channels_internal.h"

void aic_channel_mix_unitaries(aic_ucp_kraus *phi, const acb_mat_t *U,
                               const double *probs, slong m, slong d, slong prec)
{
    if (U == NULL || m < 1 || d < 1) {
        flint_printf("aic_channel_mix_unitaries: U NULL or m<1 or d<1 "
                     "(Rule 4)\n");
        flint_abort();
    }
    aic_chan_validate_probs(probs, m, AIC_CHAN_TOL, "aic_channel_mix_unitaries");
    for (slong j = 0; j < m; j++)
        aic_chan_validate_unitary(U[j], d, AIC_CHAN_TOL,
                                  "aic_channel_mix_unitaries", j, prec);

    /* K_j = sqrt(p_j) U_j^dag  (observable action sum_j p_j U_j X U_j^dag). */
    aic_ucp_kraus_init(phi, d, d, m);
    arb_t s;
    arb_init(s);
    for (slong j = 0; j < m; j++) {
        acb_mat_conjugate_transpose(phi->K[j], U[j]);   /* U_j^dag */
        arb_set_d(s, sqrt(probs[j]));
        acb_mat_scalar_mul_arb(phi->K[j], phi->K[j], s, prec);
    }
    arb_clear(s);
}

void aic_channel_group_twirl(aic_ucp_kraus *phi, const acb_mat_t *U, slong m,
                             slong d, slong prec)
{
    if (m < 1) {
        flint_printf("aic_channel_group_twirl: m<1 (Rule 4)\n");
        flint_abort();
    }
    /* Uniform weights 1/m; reuse mix_unitaries (it validates unitarity too). */
    double *probs = flint_malloc((size_t) m * sizeof(double));
    if (probs == NULL) {
        flint_printf("aic_channel_group_twirl: OOM\n");
        flint_abort();
    }
    for (slong j = 0; j < m; j++) probs[j] = 1.0 / (double) m;
    aic_channel_mix_unitaries(phi, U, probs, m, d, prec);
    flint_free(probs);
}

/* sigma_k entry setters (k = 0:I, 1:X, 2:Y, 3:Z) on a 2x2, scaled by `s`. */
static void pauli_set(acb_mat_t K, int k, double s)
{
    switch (k) {
    case 0: /* I */
        acb_set_d(acb_mat_entry(K, 0, 0), s);
        acb_set_d(acb_mat_entry(K, 1, 1), s);
        break;
    case 1: /* X = [[0,1],[1,0]] */
        acb_set_d(acb_mat_entry(K, 0, 1), s);
        acb_set_d(acb_mat_entry(K, 1, 0), s);
        break;
    case 2: /* Y = [[0,-i],[i,0]] */
        acb_set_d_d(acb_mat_entry(K, 0, 1), 0.0, -s);
        acb_set_d_d(acb_mat_entry(K, 1, 0), 0.0, s);
        break;
    case 3: /* Z = [[1,0],[0,-1]] */
        acb_set_d(acb_mat_entry(K, 0, 0), s);
        acb_set_d(acb_mat_entry(K, 1, 1), -s);
        break;
    default:
        flint_printf("pauli_set: bad index %d\n", k);
        flint_abort();
    }
}

void aic_channel_pauli(aic_ucp_kraus *phi, const double probs[4], slong prec)
{
    (void) prec;
    aic_chan_validate_probs(probs, 4, AIC_CHAN_TOL, "aic_channel_pauli");
    /* K_k = sqrt(p_k) sigma_k (Hermitian, observable = dual). Phi(X) =
     * sum_k p_k sigma_k X sigma_k; unital + trace-preserving since the sigma_k
     * are unitary and Hermitian. d = 2 (qubit). */
    aic_ucp_kraus_init(phi, 2, 2, 4);
    for (int k = 0; k < 4; k++)
        pauli_set(phi->K[k], k, sqrt(probs[k]));
}
