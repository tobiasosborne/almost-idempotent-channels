# eta0_oracle.jl — README plot 4 (Part D.2 script 4).
#
# On an EXACT conditional expectation (η = 0) every certified defect — the
# associator ε, the isomorphism defect, and both round-trip brackets — collapses to
# machine ε (exactly 0 or ~1e-36 at the arb working precision), far below the 1e-9
# reference, and the extracted block structure is EXACT: identity(2)→[2],
# block_cond_exp(4,2)→[2,2], dephasing(3)→[1,1,1]. The cleanest cross-check rung
# (test/test_factorize.jl). Grouped bar chart, log-y. Palette: Part E. ONE 640 px PNG.

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

# ---- η=0 oracle fixtures (inlined from test/fixtures.jl) ----
identity_kraus(d) = [Matrix{ComplexF64}(I, d, d)]
function block_cond_exp_kraus(d, m)
    K0 = zeros(ComplexF64, d, d); K1 = zeros(ComplexF64, d, d)
    for i in 1:m; K0[i, i] = 1; end
    for i in (m + 1):d; K1[i, i] = 1; end
    return Matrix{ComplexF64}[K0, K1]
end
function dephasing_kraus(d)
    K = Matrix{ComplexF64}[]
    for i in 1:d
        Ki = zeros(ComplexF64, d, d); Ki[i, i] = 1; push!(K, Ki)
    end
    return K
end

oracles = (("identity(2)\n→ [2]",          identity_kraus(2),       "[2]"),
           ("block_cond_exp(4,2)\n→ [2,2]",  block_cond_exp_kraus(4, 2), "[2,2]"),
           ("dephasing(3)\n→ [1,1,1]",       dephasing_kraus(3),      "[1,1,1]"))

# the four certified defects per oracle, in plot order
defect_names = ["ε (associator)", "iso defect", "‖ΔΥ−Φ‖_cb", "‖ΥΔ−1_B‖_cb"]
DISPLAY_FLOOR = 1e-40   # exact-zero / sub-floor defects render at this floor (log-y)

# collect raw values, assert the oracle invariant, build a display value with floor
raw = Matrix{Float64}(undef, length(oracles), 4)
disp = Matrix{Float64}(undef, length(oracles), 4)
blocklists = String[]
for (oi, (name, kr, blstr)) in enumerate(oracles)
    Φ = UCPMap(kr)
    A = associated_algebra(Φ; prec=128)
    v = main_isomorphism(Φ; prec=128)
    F = factorize(Φ; prec=128)
    vals = [maximum(epsilon(A)), maximum(isodefect(v)),
            maximum(delups_defect(F)), maximum(upsdel_defect(F))]
    # invariant: η=0 oracle ⇒ every defect < 1e-9, and the block list is exact
    @assert all(x -> x < 1e-9, vals) "η=0 oracle defect ≥ 1e-9 for $name: $vals"
    expected_blocks = Int[parse(Int, s) for s in split(strip(blstr, ['[', ']']), ",")]
    @assert blocks(algebra(F)) == expected_blocks "block list mismatch for $name: got $(blocks(algebra(F))), expected $expected_blocks"
    raw[oi, :] = vals
    disp[oi, :] = max.(vals, DISPLAY_FLOOR)
    push!(blocklists, blstr)
    println("$name  raw=$vals  blocks=$(blocks(algebra(F)))")
end

# ---- grouped bar chart ----
ngroups = length(oracles); nbars = 4
bar_cols = ["#7dcfff", "#9ece6a", "#bb9af7", "#e0af68"]  # palette cycle (first 4)

# flatten into (x position, height, group color) for Makie barplot with dodge
xs = Int[]; hs = Float64[]; grp = Int[]; dodge = Int[]
for oi in 1:ngroups, bi in 1:nbars
    push!(xs, oi); push!(hs, disp[oi, bi]); push!(grp, bi); push!(dodge, bi)
end

fig = Figure(size = (640, 450))
ax = Axis(fig[1, 1];
    title = "η = 0 oracle: every certified defect collapses to machine ε  (test/test_factorize.jl)",
    titlesize = 12.5,
    xlabel = "exactly-idempotent channel  (→ extracted C* block list)",
    ylabel = "max certified defect",
    yscale = log10,
    xticks = (1:ngroups, [o[1] for o in oracles]),
    xticklabelsize = 12)

ylims!(ax, DISPLAY_FLOOR / 3, 1e-6)

bp = barplot!(ax, xs, hs;
    dodge = dodge, color = [bar_cols[g] for g in grp],
    strokecolor = "#24283b", strokewidth = 0.5, gap = 0.15, dodge_gap = 0.04)

# the 1e-9 machine reference line (dashed red), the test threshold
href = hlines!(ax, [1e-9]; color = "#f7768e", linestyle = :dash, linewidth = 2)
text!(ax, 0.55, 1e-9; text = "1e-9 reference (test threshold)",
      align = (:left, :bottom), color = "#f7768e", fontsize = 12, offset = (0, 2))

# legend for the four defect types (use coloured poly elements), in the empty mid
# band well above the floored bars so it does not collide with them
deflegs = [PolyElement(color = bar_cols[i]) for i in 1:4]
axislegend(ax, deflegs, defect_names;
    position = (0.5, 0.62), framevisible = true, backgroundcolor = ("#1f2335", 0.85),
    labelcolor = "#c0caf5", patchsize = (16, 14), labelsize = 11,
    orientation = :horizontal, nbanks = 2)

# note that exact-zero bars render at the display floor (top-left, clear of bars)
text!(ax, 0.55, 1e-15; text = "all defects ε=iso=0 exactly; round-trips ≤ 1.5e-36\n(exact-0 / sub-1e-40 bars shown at the 1e-40 floor)",
      align = (:left, :top), color = "#565f89", fontsize = 11, offset = (0, 0))

out = "/home/tobias/Projects/almost-idempotent-channels/docs/assets/eta0_oracle.png"
save(out, fig; px_per_unit = 2)
println("WROTE ", out)
