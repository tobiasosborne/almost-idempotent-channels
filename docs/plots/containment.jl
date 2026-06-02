# containment.jl — README plot 1 (Part D.2 script 1).
#
# The certified eig-free bracket [lo, hi] contains the analytic idempotency defect
# η = t(1−t)·2(1−1/d²) for the depolarizing family Φ_t = (1−t)·id + t·Dep_d, over
# t ∈ [0.02, 0.30] and d ∈ {2, 3}. Containment is a test invariant
# (test/test_core.jl). Palette: Part E (Tokyo-Night-Storm). Writes ONE 640 px PNG.

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

# ---- fixture (inlined from test/fixtures.jl) ----
e(d, i) = (v = zeros(ComplexF64, d); v[i] = 1; v)
function phit_kraus(d, t)
    K = Matrix{ComplexF64}[sqrt(1 - t) * Matrix{ComplexF64}(I, d, d)]
    for i in 1:d, j in 1:d
        push!(K, sqrt(t) * (e(d, i) * e(d, j)') / sqrt(d))
    end
    return K
end

analytic_eta(d, t) = t * (1 - t) * 2 * (1 - 1 / d^2)

ts = collect(range(0.02, 0.30; length = 15))
dims = (2, 3)

# colours: band green (certified), analytic line fg, per-d marker hues
band_col = ("#9ece6a", 0.25)
line_col = "#c0caf5"
d_hue = Dict(2 => "#7dcfff", 3 => "#bb9af7")

data = Dict{Int,NamedTuple}()
for d in dims
    los = Float64[]; his = Float64[]; etas = Float64[]
    for t in ts
        Φ = UCPMap(phit_kraus(d, t))
        b = certified_defect(Φ)
        push!(los, minimum(b)); push!(his, maximum(b))
        push!(etas, analytic_eta(d, t))
        # invariant: the analytic η must lie inside the certified bracket
        @assert minimum(b) - 1e-9 ≤ analytic_eta(d, t) ≤ maximum(b) + 1e-9 "containment failed at d=$d t=$t"
    end
    data[d] = (los = los, his = his, etas = etas)
    println("d=$d: lo[1]=$(los[1])  hi[1]=$(his[1])  eta[1]=$(etas[1])  |  lo[end]=$(los[end]) hi[end]=$(his[end]) eta[end]=$(etas[end])")
end

fig = Figure(size = (640, 420))
ax = Axis(fig[1, 1];
    title = "Certified eig-free bracket contains the analytic η  (test/test_core.jl)",
    titlesize = 14,
    xlabel = "t  (mixing parameter of Φ_t = (1−t)·id + t·Dep_d)",
    ylabel = "idempotency defect  ‖Φ²−Φ‖_cb")

# shaded certified band per dimension (behind), then the analytic line on top
for d in dims
    nt = data[d]
    band!(ax, ts, nt.los, nt.his; color = band_col)
end
lineplots = []
for d in dims
    nt = data[d]
    lp = lines!(ax, ts, nt.etas; color = d_hue[d], linewidth = 2.5)
    scatter!(ax, ts, nt.etas; color = d_hue[d], markersize = 7)
    push!(lineplots, lp)
end

# legend: band entry + per-d analytic lines
band_entry = PolyElement(color = band_col)
axislegend(ax,
    [band_entry, lineplots[1], lineplots[2]],
    ["certified bracket [lo, hi]", "analytic η, d=2", "analytic η, d=3"];
    position = :lt, framevisible = false, labelcolor = "#c0caf5", patchsize = (22, 14))

out = "/home/tobias/Projects/almost-idempotent-channels/docs/assets/containment.png"
save(out, fig; px_per_unit = 2)
println("WROTE ", out)
