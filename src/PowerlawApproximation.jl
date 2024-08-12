module PowerlawApproximation

using AuditorySignalUtils
using CairoMakie
using Colors
using ColorSchemes
using Optim
using Trapz

include("profile.jl")

export pl, e, pea, pea2, fig1, fig2, loss, calc_optim

# Compute power-law function `pl` or exponential function `e`
pl(t, β) = β / (t + β)
e(t, τ) = exp(-t/τ)

# Compute sum of exponentials to approximate power-law `pea(t, τ, w)` of the form:
#   ∑ᵢ wᵢ exp(-t/τᵢ)
pea(t, τ, w=ones(length(τ))) = sum(map(x -> x[2] .* e.(t, x[1]), zip(τ, w)))
pea_components(t, τ, w=ones(length(τ))) = map(x -> x[2] .* e.(t, x[1]), zip(τ, w))

# Compute simple squared error loss computed with trapezoidal integration
function loss(t, x, y)
    trapz((t), (y .- x) .^ 2)
end

# Numerically optimize ζ parameter (SCHEME 1)
# We approximate the power-law kernel with exponential kernels with time constants τ
# and weights w where w = (1/τ)^ζ. ζ is an unknown adjustment factor that adjusts weights
# to better match the power-law kernel. Here, it is optimized numerically with Optim.jl
# and then the used τ, final w, and ζ are reported. 
function calc_optim(β=1e-2; base=10.0, start=0, step=1/exp(1), stop=3, fs=10e3)
    # Determine time constants
    τ = β .* base .^ collect(start:step:stop)

    # Set initial parameters and bounds
    ζ = 0.2
    w = 1 ./ τ

    # Synthesize time vector and pl kernel
    t = timevec(10.0, fs)
    kernel = pl.(t, β)

    # Define loss function
    f = x -> loss(t, kernel, β .* sum(pea_components(t, τ, w .^ x[1]))) 

    # Compute minimization
    ζ̂ = Optim.minimizer(optimize(f, [ζ]))
    
    # Return results (for convenience, return τ, w, and ζ)
    τ, w .^ ζ̂, ζ̂
end

# Figure 1
# Show target power-law kernel and approximation, as well as each component contributing
# to the fit, on a log-log scale.
function fig1(
    τ,
    w;
    β=1e-2, 
    fs=10e3, 
    xscale=log10, 
    yscale=log10,
    xlabel="Time (s)",
    ylabel="Amplitude (a.u.)",
    size=(1000, 700),
    dur=1e1,
    ylims=(1e-5, 5e0),
    xlims=(1/fs / 10, dur * 10),
    fig=Figure(; size=size), 
    ax=Axis(
        fig[1, 1]; 
        xscale=xscale,
        yscale=yscale,
        xlabel=xlabel,
        ylabel=ylabel,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
    ),
    colorscheme=:glasgow,
)
    # Create time vector & compute PEA
    t = timevec(dur, fs)
    g_comp = β .* pea_components(t, τ, w)
    g = β .* pea(t, τ, w)
    f = pl.(t, β)

    # Select colors for each τᵢ
    colors = get(colorschemes[colorscheme], LinRange(0.0, 1.0, length(τ)))

    # Start figure with reference lines:
    #   1) horizontal gridline at 0.5
    hlines!(ax, [0.5]; linestyle=:dash, color=:gray)

    # Plot results:
    #   1) Approximation components in colored lines, matching to legend
    #   2) Composite approximation in pink
    #   3) Target in red
    map(zip(g_comp, τ, colors)) do (gᵢ, τᵢ, colorᵢ)
        lines!(ax, t[2:end], gᵢ[2:end]; color=colorᵢ, label=string(round(τᵢ*1e3)) * " ms")
    end
    lines!(ax, t[2:end], g[2:end]; color=:pink, linestyle=:dash)
    lines!(ax, t[2:end], f[2:end]; color=:red)

    # Adjust limits
    ylims!(ax, ylims)
    xlims!(ax, xlims)
    
    # Add legend
    axislegend(ax)
    fig
end

end # module PowerlawApproximation
