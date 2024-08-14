# genfigs.jl
# 
# Generate figures for Guest and Carney (2024), JASA
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
fig = fig1([1e-2], [1/β]; β=β, size=(400, 300))
saveplot("fig1a.png", fig)
fig = fig1([1e-2], [1/β]; β=β, size=(400, 300), xscale=identity, yscale=identity, ylims=(0.0, 1.0), xlims=(0.0, 5β))
saveplot("fig1b.png", fig)

# Figure 2?
fig = fig1()
saveplot("fig1.png", fig)
