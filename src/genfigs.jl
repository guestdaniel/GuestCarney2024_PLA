# genfigs.jl
# 
# Generate figures for Guest and Carney (2024), JASA
#
# A NOTE ON CODE ORGANIZATION:
# The "theory" and numerical optimization code is all written in Julia and conveniently
# provided by this package. The "testing" code where the approximation scheme is tested
# in practice is written in MATLAB because the target wrapper system is Mex. That code
# is included here in the folder `matlab`. The underlying model code, of course, is written
# in C, and is included in `src\models`.
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

# Figure 1
β = 1e-2
sz = (250, 200)
fig = fig1([1e-2], [1/β]; β=β, size=sz, xlims=(5e-4, 1.1e1), ylims=(5e-4, 5.0))
vlines!([β]; color=:gray)
text!([β*1.1], [1.0]; text="β=$β s", color=:gray)
saveplot("fig1a.png", fig)
fig = fig1([1e-2], [1/β]; β=β, size=sz, xscale=identity, yscale=identity, ylims=(0.0, 1.0), xlims=(0.0, 5β))
vlines!([β]; color=:gray)
text!([β + 1e-3], [0.85]; text="β=$β s", color=:gray)
saveplot("fig1b.png", fig)

# Figure 2?
β = 1e-2
τ = β .* 10 .^ (0.0:0.5:3)
w = 1 ./ τ
fig = fig1(τ, w; size=(sz[1]*1.7, sz[2]*1.7), xlims=(1e-4, 1e2), dur=50.0)
saveplot("fig1c.png", fig)
