/* aic_idemp_wedderburn_assemble.c — the W_j assembly of the Artin-Wedderburn
 * decomposition (bead aic-ynu, I1): group eigen-clusters into Wedderburn blocks
 * by W-connectivity, then build each block-diagonalizing isometry W_j by aligning
 * a single E_j onb across the dim_L_j reference vectors. Double path; the output
 * W_j is certified in arb (dim identities + W-unitarity + the w-reconstruction
 * correctness spec) by aic_wedd_certify_W BEFORE this routine returns, so a
 * tensor-MISALIGNED W is rejected in PRODUCTION. See aic_idemp_wedderburn.c for the
 * full narrative (.tex:259-272, :2106). The small linear-algebra helpers +
 * clustering + the certification live in aic_idemp_wedderburn_ops.c.
 *
 * THE BLOCKS. Q_a, Q_a' are in the SAME Wedderburn block iff the algebra connects
 * them: max_k ||Q_a w(B_k) Q_a'||_op > tol (every w(B_k) is block-diagonal across
 * j, so cross-block products vanish). Blocks = connected components (union-find).
 *
 * RE-DRAW (the aic-ynu BLOCKER fix). A degenerate H_W draw (a cluster spanning two
 * would-be blocks -> disagreeing in-block cluster multiplicities, or a vanishing
 * connecting isometry) is NOT an out-of-scope input: it is a measure-zero accident
 * of THIS random H_W. We return nonzero (out reset to empty) so the driver re-draws
 * a fresh H_W; only the last draw in the budget fails loud (Rule 4).
 *
 * W_j ASSEMBLY (the tensor alignment — load-bearing). Block j has reference
 * cluster a_0 (onb {f_c}, c in E_j) and dim_L_j clusters a_p. The connecting map
 * T_p = sum_k Q_{a_p} w(B_k) Q_{a_0} sends v_0 (x) E_j -> v_p (x) E_j as
 * (v_p v_0^dag) (x) C; applying it to f_c and projecting onto cluster a_p's span,
 * then normalizing, yields an E_j onb at v_p ALIGNED (same index c) to v_0's. The
 * W_j row t = a*dim_E + c (L-MAJOR, matching aic_mat_kronecker A_j (x) 1_{E_j}) is
 * the bra f_{a,c}^dag, i.e. W_j[t, :] = conj(f_{a,c}).
 */
#include <assert.h>
#include <complex.h>
#include <math.h>

#include <flint/acb_mat.h>

#include "aic/aic_latd.h"
#include "aic_idemp_wedderburn_internal.h"

/* max_k ||Q_a w(B_k) Q_b||_op, the W-connectivity of two clusters (lives here, with
 * the union-find that consumes it). Q_a Wk Q_b = F_a G_k F_b^dag with
 * G_k[ca,cb] = <f_{a,ca}, Wk f_{b,cb}>; ||Q_a Wk Q_b||_op = ||G_k||_op (F_a, F_b
 * isometries). Zero across distinct Wedderburn blocks. */
double aic_wedd_block_link(const double _Complex *evecs, slong sa, slong ka,
                           slong sb, slong kb, double _Complex *const *Wk,
                           slong dA, slong m)
{
    double best = 0.0;
    double _Complex *col = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *Wf = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *fa = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *G = flint_malloc((size_t) (ka * kb) * sizeof(double _Complex));
    for (slong k = 0; k < dA; k++) {
        for (slong cb = 0; cb < kb; cb++) {
            aic_wedd_getcol(col, evecs, sb + cb, m);
            aic_wedd_matvec(Wf, Wk[k], col, m);
            for (slong ca = 0; ca < ka; ca++) {
                aic_wedd_getcol(fa, evecs, sa + ca, m);
                G[ca * kb + cb] = aic_wedd_cdot(fa, Wf, m);
            }
        }
        double s = aic_latd_opnorm(G, ka, kb);
        if (s > best) best = s;
    }
    flint_free(G);
    flint_free(fa);
    flint_free(Wf);
    flint_free(col);
    return best;
}

/* write the bra f^dag into row t of W (W[t,i] = conj(f_i)). */
static void set_row_bra(acb_mat_t W, slong t, const double _Complex *f, slong m)
{
    for (slong i = 0; i < m; i++)
        acb_set_d_d(acb_mat_entry(W, t, i), creal(conj(f[i])), cimag(conj(f[i])));
}

/* Free any W_j already allocated and reset `out` to empty. Used when a degenerate
 * H_W draw must be abandoned for a re-draw (see the RE-DRAW PROTOCOL). */
static void wedd_reset_out(aic_idemp_wedderburn *out, slong built)
{
    if (out->W_j != NULL) {
        for (slong j = 0; j < built; j++) acb_mat_clear(out->W_j[j]);
        flint_free(out->W_j);
        out->W_j = NULL;
    }
    flint_free(out->dim_L); out->dim_L = NULL;
    flint_free(out->dim_E); out->dim_E = NULL;
    out->num_blocks = 0;
}

int aic_wedd_assemble(aic_idemp_wedderburn *out, double _Complex *const *Wk,
                      slong dA, const double _Complex *evecs,
                      const slong *cl_start, const slong *cl_size, slong n_cl,
                      slong m, slong prec, int fail_loud)
{
    /* connectivity: union-find over clusters via block_link > tol. At eta=0 the
     * in-block vs cross-block margin is sharp and O(1): same-block clusters connect
     * with ||Q_a w(B_k) Q_b||_op = O(1) (the algebra acts transitively within a
     * Wedderburn block), while DISTINCT blocks have it EXACTLY 0 (every w(B_k) is
     * block-diagonal), floored only by the double-path zheev backward error ~1e-14.
     * link_tol=1e-7 sits in that O(1)-to-1e-14 gulf, safe with ~7 orders of slack;
     * not prec-scaled because the margin is set by the exact eta=0 block structure,
     * not by `prec`. (A perturbed eta>0 input is out of I1 scope.) */
    slong *comp = flint_malloc((size_t) n_cl * sizeof(slong));
    for (slong a = 0; a < n_cl; a++) comp[a] = a;
    const double link_tol = 1e-7;
    for (slong a = 0; a < n_cl; a++)
        for (slong b = a + 1; b < n_cl; b++) {
            double l = aic_wedd_block_link(evecs, cl_start[a], cl_size[a],
                                           cl_start[b], cl_size[b], Wk, dA, m);
            if (l > link_tol) {
                slong ra = a, rb = b;
                while (comp[ra] != ra) ra = comp[ra];
                while (comp[rb] != rb) rb = comp[rb];
                if (ra < rb) comp[rb] = ra; else comp[ra] = rb;
            }
        }
    for (slong a = 0; a < n_cl; a++) {
        slong r = a; while (comp[r] != r) r = comp[r];
        comp[a] = r;
    }

    /* distinct roots -> blocks, in order of first appearance. */
    slong *root_of = flint_malloc((size_t) n_cl * sizeof(slong));
    slong nb = 0;
    for (slong a = 0; a < n_cl; a++) {
        slong found = -1;
        for (slong j = 0; j < nb; j++) if (root_of[j] == comp[a]) { found = j; break; }
        if (found < 0) { root_of[nb] = comp[a]; nb++; }
    }

    out->num_blocks = nb;
    out->dim_L = flint_malloc((size_t) nb * sizeof(slong));
    out->dim_E = flint_malloc((size_t) nb * sizeof(slong));
    out->W_j = flint_malloc((size_t) nb * sizeof(acb_mat_t));

    double _Complex *f0 = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *Wf = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *acc = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *proj = flint_malloc((size_t) m * sizeof(double _Complex));
    double _Complex *fa = flint_malloc((size_t) m * sizeof(double _Complex));

    int redraw = 0;          /* set if a degenerate draw is detected (retryable) */
    slong built = 0;         /* number of W_j[] already acb_mat_init'd           */
    for (slong j = 0; j < nb && !redraw; j++) {
        /* clusters in this block, ascending eigenvalue order (= array order). */
        slong *cls = flint_malloc((size_t) n_cl * sizeof(slong));
        slong nL = 0;
        for (slong a = 0; a < n_cl; a++)
            if (comp[a] == root_of[j]) cls[nL++] = a;
        slong dE = cl_size[cls[0]];
        /* A genuine block (M_d simple or M_d (+) M_d) has every in-block cluster of
         * the SAME multiplicity dim_E. Disagreeing cluster sizes within a "block"
         * mean a degenerate H_W draw landed two distinct blocks in one cluster
         * (the aic-ynu BLOCKER) — re-drawable, NOT an out-of-scope input. */
        for (slong p = 0; p < nL; p++)
            if (cl_size[cls[p]] != dE) { redraw = 1; break; }
        if (redraw) { flint_free(cls); break; }
        out->dim_L[j] = nL;
        out->dim_E[j] = dE;
        acb_mat_init(out->W_j[j], nL * dE, m);
        built = j + 1;

        for (slong c = 0; c < dE && !redraw; c++) {
            aic_wedd_getcol(f0, evecs, cl_start[cls[0]] + c, m); /* reference f_c */
            set_row_bra(out->W_j[j], 0 * dE + c, f0, m);          /* p = 0 row */
            for (slong p = 1; p < nL; p++) {
                slong a = cls[p];
                /* T_p f_c = sum_k Wk f_c, projected onto cluster a's span. */
                for (slong i = 0; i < m; i++) acc[i] = 0;
                for (slong k = 0; k < dA; k++) {
                    aic_wedd_matvec(Wf, Wk[k], f0, m);
                    for (slong i = 0; i < m; i++) acc[i] += Wf[i];
                }
                for (slong i = 0; i < m; i++) proj[i] = 0;
                for (slong cc = 0; cc < dE; cc++) {
                    aic_wedd_getcol(fa, evecs, cl_start[a] + cc, m);
                    double _Complex coeff = aic_wedd_cdot(fa, acc, m);
                    for (slong i = 0; i < m; i++) proj[i] += coeff * fa[i];
                }
                double nrm = sqrt(creal(aic_wedd_cdot(proj, proj, m)));
                /* A vanishing connecting isometry: sum_k w(B_k) failed to carry the
                 * reference E_j onto cluster a — again a degenerate draw (the
                 * generator combination annihilated the link). Re-drawable.
                 * Threshold 1e-9: a genuine connecting isometry has unit-O(1) norm
                 * (it maps an onb vector to an onb vector at eta=0), so the gap from
                 * a true link to a numerically-zero one is O(1) vs ~1e-14; 1e-9 is
                 * deep in that gulf, matching link_tol's margin. */
                if (nrm < 1e-9) { redraw = 1; break; }
                for (slong i = 0; i < m; i++) proj[i] /= nrm;
                set_row_bra(out->W_j[j], p * dE + c, proj, m);
            }
        }
        flint_free(cls);
    }

    flint_free(fa);
    flint_free(proj);
    flint_free(acc);
    flint_free(Wf);
    flint_free(f0);
    flint_free(root_of);
    flint_free(comp);

    if (redraw) {
        wedd_reset_out(out, built);
        if (fail_loud) {
            flint_printf("aic_idemp_wedderburn: H_W draw still degenerate after the "
                         "full re-draw budget (a cluster spans two blocks, or a "
                         "connecting isometry vanished); genuinely-unresolvable "
                         "spectrum -> fail loud (Rule 4)\n");
            flint_abort();
        }
        return 1;        /* caller re-draws a fresh H_W */
    }

    /* Production certification (arb): dim identities + W-unitarity + the
     * w-reconstruction correctness spec (rejects a tensor-MISALIGNED W, FIX 2).
     * Fails loud on any violation; a clean return here means a certified W. */
    aic_wedd_certify_W(out, Wk, dA, prec);
    return 0;
}
