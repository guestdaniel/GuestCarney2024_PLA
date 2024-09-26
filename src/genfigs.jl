# genfigs.jl
# 
# Generate figures for Guest and Carney (2024), JASA
#
# A NOTE ON CODE ORGANIZATION:
# The "theory" and numerical optimization code is all written in Julia and conveniently
# provided by this package. The "testing" code where the approximation scheme is tested
# in practice is written in MATLAB because the target wrapper system is Mex. That code
# is included here in the folder `src\matlab`. The underlying model code, of course, is written
# in C, and is included in `src\c`.
#
# A NOTE ON FIGURE SIZES:
# Plot sizes are initially specified in pixels, but implicitly we assume that because we 
# will render at 300 DPI, that plot sizes are specified in inches*100... i.e., if we take 
# plot size (600, 400)px and do px_per_unit=3 we achieve (1800, 1200)px and or (6, 4)in at 
# 300 DPI. This is mostly handled by the `saveplot` function below.

# Load packages
using CairoMakie
using PowerlawApproximation

# Next, define a convenience function for saving plots
saveplot(name, fig) = save(projectdir("figs", name), fig; px_per_unit=3) 

# Next, set the theme used for the figures 
set_theme!()

# Next, set some settings that get used throughout
sz = (250, 200)  # basic figure size, 2.5" x 2"

# Figure 1
# Figure depicting the basic principles of approximating power-law kernel with a weighted
# sum of exponential kernels
β = 1e-2

fig = fig1([1e-2], [1/β]; β=β, size=sz, xlims=(5e-4, 1.1e1), ylims=(5e-4, 5.0))
vlines!([β]; color=:gray)
text!([β*1.1], [1.0]; text="β=$β s", color=:gray)
saveplot("fig1a.png", fig)

fig = fig1([1e-2], [1/β]; β=β, size=sz, xscale=identity, yscale=identity, ylims=(0.0, 1.0), xlims=(0.0, 5β))
vlines!([β]; color=:gray)
text!([β + 1e-3], [0.85]; text="β=$β s", color=:gray)
saveplot("fig1b.png", fig)

τ = β .* 10 .^ (0.0:0.5:3)
w = 1 ./ τ
fig = fig1(τ, w; size=(sz[1]*1.7, sz[2]*1.7), xlims=(1e-4, 1e2), dur=50.0)
saveplot("fig1c.png", fig)

# Fig 2
# Comparing optimized and heuristic weights for approximation
τ = β .* 10 .^ (0:1/exp(1):5)
w = 1 ./ (τ .+ β)
c = 2β * sum(w .* exp.(-β ./ τ))
w = w ./ c
fig = fig1(τ, w; size=(sz[1]*1.7, sz[2]*1.7), xlims=(1e-4, 1e2), dur=50.0)
saveplot("fig2a.png", fig)

fig = fig1(calc_optim_s2(; start=-1, stop=5, dur=50.0)...; size=(sz[1]*1.7, sz[2]*1.7), xlims=(1e-4, 1e2), dur=50.0)
saveplot("fig2b.png", fig)

# Fig 3
# Actual weights/taus used for β values found in ZBC2014 model
fig = fig1(
    calc_heuristic_s1(5e-4; start=0, stop=5)...; 
    β=5e-4,
    size=(sz[1], sz[2]), 
)
saveplot("fig3a_heuristic_slow.png", fig)

fig = fig1(
    calc_heuristic_s1(1e-1; start=0, stop=5)...; 
    β=1e-1,
    size=(sz[1], sz[2]), 
)
saveplot("fig3a_heuristic_fast.png", fig)

fig = fig1(
    calc_optim_s2(5e-4; start=0, stop=5, dur=5e-4*1e5)...;
    β=5e-4,
    size=(sz[1], sz[2]), 
)
saveplot("fig3a_optimized_slow.png", fig)

fig = fig1(
    calc_optim_s2(1e-1; start=0, stop=5, dur=1e-1*1e5)...; 
    β=1e-1,
    size=(sz[1], sz[2]), 
)
saveplot("fig3a_optimized_fast.png", fig)