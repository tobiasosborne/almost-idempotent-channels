# test_factorize.jl — the SOLVER-FREE factorize round-trip suite (bead aic-exa.9
# [T], design §9 step 2 / Appendix B1/B4). The headline verb: Φ ≈ Δ Υ through a
# genuine C* algebra B, with the two certified cb-norm round-trip brackets the
# rigorous certificate (NOT iscptp at machine tolerance — Appendix B4).
#
# THE η=0 IDENTITY ORACLE IS STRUCTURALLY BLIND TO STAR/DUAL BUGS (§C2/§C13): an
# exactly-idempotent Φ has Φ̃(XY) = XY on its image, so a "plain XY" star bug and a
# Dec/Enc dual SWAP both read 0 there. So the OBLIQUE η>0 fixture (bce_conj /
# mixconj) is MANDATORY — it is the only fixture that exercises the Choi–Effros star
# and the dual-direction binding. We assert:
#   (i)  η=0 oracle: delups/upsdel straddle 0 at machine-ε; B-blocks match the
#        aic_idemp structure (identity→[d], block_cond_exp(4,2)→[2,2]).
#   (ii) η>0 oblique: delups/upsdel are O(η) (not 0, not O(1)); encode/decode are
#        CPMaps with the B1 dims ASSERTED (so a swap is caught); the MULTI-BLOCK
#        oblique decode is only O(η)-TP (iscptp@1e-9 false, @1e-4 true) — the §C14
#        PSD-cone clip; the round-trip bracket is the real certificate.
#
# Assumes `using AlmostIdempotentChannels, Test, LinearAlgebra` + fixtures.jl.

@testset "η=0 ORACLE: identity / exact-idempotent factorize (machine-ε straddle)" begin
    for (name, kr, wantblocks) in (("identity(2)", identity_kraus(2), [2]),
                                   ("block_cond_exp(4,2)", block_cond_exp_kraus(4, 2), [2, 2]),
                                   ("dephasing(3)", dephasing_kraus(3), [1, 1, 1]))
        Φ = UCPMap(kr)
        F = factorize(Φ)
        # B-blocks match the aic_idemp (Wedderburn) structure of the exact idempotent.
        @test sort(blocks(algebra(F))) == sort(wantblocks)
        # delups = ‖ΔΥ − Φ‖_cb and upsdel = ‖ΥΔ − 1_B‖_cb straddle 0 at machine-ε.
        @test minimum(delups_defect(F)) <= 0.0 <= maximum(delups_defect(F))
        @test minimum(upsdel_defect(F)) <= 0.0 <= maximum(upsdel_defect(F))
        @test maximum(delups_defect(F)) < 1e-9
        @test maximum(upsdel_defect(F)) < 1e-9
        @test F.eta_proxy < 1e-9              # η=0 sentinel
        # η=0 decode is EXACTLY trace-preserving (no PSD-cone clip).
        @test iscptp(decode(F); atol=1e-9)
        @test iscptp(encode(F); atol=1e-9)
        @info "η=0 factorize oracle" name blocks=blocks(algebra(F)) delups_hi=maximum(delups_defect(F)) upsdel_hi=maximum(upsdel_defect(F))
    end
end

@testset "η>0 OBLIQUE factorize: dims + O(η) round-trip + O(η)-TP decode (B1/B4)" begin
    # bce_conj(4,2,t): block_cond_exp(4,2) conjugate-mixed → a MULTI-BLOCK (B = M_2 ⊕
    # M_2) oblique η>0 channel. The genuine multi-block oblique structure triggers the
    # §C14 Choi→Kraus PSD-cone clip, so the decode is only O(η)-TP. This is the
    # MANDATORY §C2/§C13 fixture — the η=0 oracle is blind to the star/dual.
    for t in (0.02, 0.05)
        Φ = UCPMap(bce_conj_kraus(4, 2, t); check=true)
        @test isunital(Φ)
        F = factorize(Φ)
        η = F.eta_proxy
        @test η > 1e-3                       # genuinely η>0
        B = algebra(F)
        @test sort(blocks(B)) == [2, 2]       # multi-block oblique B = M_2 ⊕ M_2
        N = ambient(associated_algebra(Φ))    # = 4
        # B1 DIMS — assert so a Dec/Enc swap is caught: encode = Enc = Υ* : B→B(H)
        # (dim_in = n_B, dim_out = N); decode = Dec = Δ* : B(H)→B (dim_in = N, dim_out = n_B).
        en = encode(F)
        de = decode(F)
        @test en.dim_in == n_B(B)
        @test en.dim_out == N
        @test de.dim_in == N
        @test de.dim_out == n_B(B)
        # the encode/decode are DUALS: encode dims = reverse of decode dims (a swap
        # would make them equal, NOT reversed). (Here n_B = 2+2 = 4 = N, so this
        # fixture's CPMaps are square; the rectangular case is the single-block
        # mixconj testset below where n_B = 2 ≠ N = 4.)
        @test (en.dim_in, en.dim_out) == (de.dim_out, de.dim_in)
        # delups/upsdel are O(η): NOT 0 (would mean a degenerate factorization) and
        # NOT O(1) (would mean a broken one). The bracket is the rigorous certificate.
        @test maximum(delups_defect(F)) > 1e-6
        @test maximum(upsdel_defect(F)) > 1e-6
        @test maximum(delups_defect(F)) < 50 * η     # O(η), generous constant (§D2)
        @test maximum(upsdel_defect(F)) < 50 * η
        # Appendix B4: decode is only O(η)-TP (the §C14 PSD-cone clip). iscptp@1e-9 is
        # FALSE; iscptp@1e-4 is TRUE; the cb-norm round-trip bracket (upsdel/delups),
        # NOT iscptp at machine tol, is the rigorous certificate. encode IS exactly TP.
        @test !iscptp(de; atol=1e-9)
        @test iscptp(de; atol=1e-4)
        @test iscptp(en; atol=1e-9)
        # cross-check the measured TP defect lies in (1e-9, 1e-4) — the O(η) clip.
        tpd = _tp_defect(de)
        @test 1e-9 < tpd < 1e-4
        @info "η>0 oblique factorize (bce_conj 4,2)" t eta_proxy=η blocks=blocks(B) delups_hi=maximum(delups_defect(F)) upsdel_hi=maximum(upsdel_defect(F)) decode_TP_defect=tpd
    end
end

@testset "η>0 OBLIQUE factorize: single-block mixconj round-trip + dims" begin
    # mixconj(4,2,0.02): compress_idemp(4,2) conjugate-mixed → a SINGLE-BLOCK (B =
    # M_2) oblique η>0 channel in factorize's small-η domain (ρ ≈ 0.04). Exercises
    # the star + dual binding on a different oblique geometry; the single-block PSD
    # cone is empty so its decode is TP to machine precision.
    Φ = UCPMap(mixconj_kraus(4, 2, 0.02); check=true)
    F = factorize(Φ)
    η = F.eta_proxy
    @test η > 1e-3
    B = algebra(F)
    @test sort(blocks(B)) == [2]
    en = encode(F)
    de = decode(F)
    @test (en.dim_in, en.dim_out) == (n_B(B), 4)      # B1 dims
    @test (de.dim_in, de.dim_out) == (4, n_B(B))
    @test maximum(delups_defect(F)) > 1e-6
    @test maximum(upsdel_defect(F)) > 1e-6
    # O(η) on BOTH round-trip directions (symmetric with the bce_conj block; the
    # earlier single-sided upsdel-only assertion was the NIT-4 asymmetry).
    @test maximum(delups_defect(F)) < 50 * η
    @test maximum(upsdel_defect(F)) < 50 * η
    @info "η>0 oblique factorize (mixconj 4,2,0.02)" eta_proxy=η blocks=blocks(B) delups_hi=maximum(delups_defect(F)) upsdel_hi=maximum(upsdel_defect(F))
end

@testset "§D2 DIM-INDEPENDENCE canary: upsdel/η, delups/η bounded + no trend in n" begin
    # The paper's constant is UNIVERSAL and dimension-INDEPENDENT (th_main, .tex:460;
    # the naive Haar/cohomology route fails precisely because its error ∝ n, .tex:484).
    # The C side has t_dim_indep_mixconj (test_shim_factorize_artifacts.c); this is the
    # mandatory Julia analogue the design §D2 calls for (NIT 3, aic-exa.9 review): sweep
    # the SAME mixconj family at t=0.02 across n=4,5,6 and assert the round-trip ratios
    # c = ‖ΥΔ−1_B‖_cb/η and ‖ΔΥ−Φ‖_cb/η stay BOUNDED and show NO upward trend in n. A
    # constant that grew with n (the failure mode the paper warns about) would go RED.
    fx = ((4, 2), (5, 3), (6, 3))
    c_ups = Float64[]
    c_del = Float64[]
    rows = NamedTuple[]
    for (d, m) in fx
        Φ = UCPMap(mixconj_kraus(d, m, 0.02); check=true)
        F = factorize(Φ)
        η = F.eta_proxy
        @test η > 1e-3
        cu = maximum(upsdel_defect(F)) / η
        cd = maximum(delups_defect(F)) / η
        @test 0.0 < cu < 20.0          # bounded (the C canary's < 20 gate)
        @test 0.0 < cd < 20.0
        push!(c_ups, cu)
        push!(c_del, cd)
        push!(rows, (n=d, c_ups=cu, c_del=cd, blocks=blocks(algebra(F))))
    end
    # NO upward trend: the largest-n ratio is not materially above the smallest-n one
    # (a per-step blow-up factor < 1.5 is generous; measured ≈ 1.0 — the ratios are
    # essentially FLAT in n, which is the dim-independence the paper claims).
    @test maximum(c_ups) <= 1.5 * minimum(c_ups)
    @test maximum(c_del) <= 1.5 * minimum(c_del)
    @info "§D2 dim-independence canary (mixconj t=0.02, n=4,5,6)" rows c_ups_range=(minimum(c_ups), maximum(c_ups)) c_del_range=(minimum(c_del), maximum(c_del))
end
