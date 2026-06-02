# roundtrip.jl — README plot 3 (Part D.2 script 3).
#
# The certified channel-factorization round-trip defects ‖ΔΥ−Φ‖_cb (delups) and
# ‖ΥΔ−1_B‖_cb (upsdel) scale LINEARLY in η for oblique multi-block channels, well
# below the generous 50·η envelope the tests assert (test/test_factorize.jl).
# Fixtures: bce_conj(4,2,t) (multi-block B = M_2 ⊕ M_2) and mixconj(4,2,t)
# (single-block B = M_2), over t ∈ {0.01,...,0.07} inside factorize's domain
# (ρ < 0.10). Log–log axes; slope ≈ 1 is the message. Palette: Part E. ONE 640 px PNG.

using CairoMakie
using AlmostIdempotentChannels
using LinearAlgebra

# ---- Part E.2 unified theme (verbatim) ----
set_theme!(Theme(
    backgroundcolor = "#24283b",
    textcolor = "#c0caf5",
    Axis = (backgroundcolor="#1f2335", xgridcolor=("#565f89",0.4),
            ygridcolor=("#565f89",0.4), xtickcolor="#565f89", ytickcolor="#565f89",
            leftspinecolor="#565f89", bottomspinecolor="#565f89",
            topspinevisible=false, rightspinevisible=false),
    palette = (color = ["#7dcfff","#9ece6a","#bb9af7","#e0af68","#f7768e","#7aa2f7"],),
    fontsize = 15,
))

# ---- fixtures (inlined from test/fixtures.jl) ----
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end
function compress_idemp_kraus(d, m)
    K = Matrix{ComplexF64}[zeros(ComplexF64, d, d)]
    for i in 1:m; K[1][i, i] = 1; end
    for b in 0:(d - m - 1)
        Kb = zeros(ComplexF64, d, d); Kb[1, m + b + 1] = 1; push!(K, Kb)
    end
    return K
end
function fixed_unitary(nn)
    H = zeros(ComplexF64, nn, nn)
    for p in 1:(nn - 1); H[p, p + 1] = 0.4; H[p + 1, p] = 0.4; end
    for p in 1:(nn - 2); H[p, p + 2] = 0.3im; H[p + 2, p] = -0.3im; end
    for p in 1:nn; H[p, p] = 0.1 * p; end
    return exp(im * Matrix(Hermitian(H)))
end
function _mix_conj_of(base, d, t)
    V = fixed_unitary(d); Vd = V'
    K = Matrix{ComplexF64}[]
    for Ka in base; push!(K, sqrt(1 - t) * Ka); end
    for Ka in base; push!(K, sqrt(t) * (Vd * Ka * V)); end
    return K
end
mixconj_kraus(d, m, t) = _mix_conj_of(compress_idemp_kraus(d, m), d, t)
bce_conj_kraus(d, m, t) = _mix_conj_of(block_cond_exp_kraus(d, m), d, t)

ts = (0.01, 0.02, 0.03, 0.05, 0.07)
families = (("bce_conj (B = M_2 ⊕ M_2)", bce_conj_kraus),
            ("mixconj (B = M_2)",        mixconj_kraus))

# collect (η, delups, upsdel) per family
results = Dict{String,NamedTuple}()
for (name, kf) in families
    ηs = Float64[]; dus = Float64[]; uds = Float64[]
    for t in ts
        Φ = UCPMap(kf(4, 2, t))
        F = factorize(Φ; prec=128)
        η = F.eta_proxy
        du = maximum(delups_defect(F)); ud = maximum(upsdel_defect(F))
        push!(ηs, η); push!(dus, du); push!(uds, ud)
        # invariant: round-trip defects sit below the 50·η envelope the tests assert
        @assert du ≤ 50 * η && ud ≤ 50 * η "50η envelope breached: $name t=$t"
    end
    results[name] = (η = ηs, delups = dus, upsdel = uds)
    println("$name: η=$(round.(ηs; sigdigits=3))  delups=$(round.(dus; sigdigits=3))  upsdel=$(round.(uds; sigdigits=3))")
end

fig = Figure(size = (640, 440))
ax = Axis(fig[1, 1];
    title = "Certified factorization round-trip is O(η)  (test/test_factorize.jl)",
    titlesize = 13.5,
    xlabel = "η  (eta_proxy = ‖Φ²−Φ‖ scale)",
    ylabel = "round-trip defect (cb-norm)",
    xscale = log10, yscale = log10)

# colour map: delups in teal, upsdel in purple; bce_conj solid circles, mixconj squares
col_delups = "#7dcfff"
col_upsdel = "#bb9af7"
markers = Dict("bce_conj (B = M_2 ⊕ M_2)" => :circle, "mixconj (B = M_2)" => :rect)

handles = []; labels = String[]
for (name, _) in families
    r = results[name]; mk = markers[name]
    h1 = scatterlines!(ax, r.η, r.delups; color = col_delups, marker = mk, markersize = 11, linewidth = 2)
    h2 = scatterlines!(ax, r.η, r.upsdel; color = col_upsdel, marker = mk, markersize = 11, linewidth = 2)
    push!(handles, h1); push!(labels, "‖ΔΥ−Φ‖_cb  ($name)")
    push!(handles, h2); push!(labels, "‖ΥΔ−1_B‖_cb  ($name)")
end

# the 50·η envelope (dashed red), the test's asserted bound
allη = vcat([results[name].η for (name, _) in families]...)
ηlo, ηhi = minimum(allη), maximum(allη)
xenv = 10 .^ range(log10(ηlo) - 0.05, log10(ηhi) + 0.05; length = 50)
henv = lines!(ax, xenv, 50 .* xenv; color = "#f7768e", linestyle = :dash, linewidth = 2)
push!(handles, henv); push!(labels, "50·η envelope (test bound)")

# reference slope-1 guide (faint) to show linear-in-η
hslope = lines!(ax, xenv, 5 .* xenv; color = ("#565f89", 0.8), linestyle = :dot, linewidth = 1.5)
push!(handles, hslope); push!(labels, "slope 1 (linear in η)")

axislegend(ax, handles, labels;
    position = :lt, framevisible = false, labelcolor = "#c0caf5",
    patchsize = (20, 12), labelsize = 11, rowgap = 1)

out = "/home/tobias/Projects/almost-idempotent-channels/docs/assets/roundtrip.png"
save(out, fig; px_per_unit = 2)
println("WROTE ", out)
