module PowerlawApproximation

using AuditorySignalUtils
using CairoMakie
using Colors
using ColorSchemes
using Optim
using Trapz
using FFTW
using DSP
using Distributed

include("proofs.jl")

# Export basic kernel/approximator functions
export pl, e, pea, pea_components
# Export figure generation functions
export fig1, fig2, fig3, fig4, fig5, fig6, fig7
# Export weight optimization functions
export loss, calc_optim_s1, calc_optim_s2
# Export heuristic weight functions
export calc_heuristic_s1
# Export Distributed setup macro
export @parallel

# Convenient macro to set up parallel compute
macro parallel(n=4)
    quote
        using Distributed
        if nprocs() == nworkers() == 1
            addprocs($n)
        elseif nworkers() < $n
            addprocs($n - nworkers())
        end
        @everywhere using Pkg
        @everywhere Pkg.activate(".")
        @everywhere using PowerlawApproximation
        @info "Parallel pool established with $(nworkers()) workers, PowerlawApproximation.jl loaded!"
    end
end

# Compute power-law function `pl` or exponential function `e`
pl(t, β) = β / (t + β)
e(t, τ) = exp(-t/τ)

# Compute sum of exponentials to approximate power-law `pea(t, τ, w)` of the form:
#   ∑ᵢ wᵢ exp(-t/τᵢ)
pea(t, τ, w=ones(length(τ))) = sum(map(x -> x[2] .* e.(t, x[1]), zip(τ, w)))
pea_components(t, τ, w=ones(length(τ))) = map(x -> x[2] .* e.(t, x[1]), zip(τ, w))

# Compute SSE loss
function loss(t, x, y; scale=identity)
    if scale == identity
        sum((y .- x).^2)
    else
        trapz((log.(t[2:end])), (log.(y[2:end]) .- log.(x[2:end])).^2)
    end
end

# calc_optim_s1 (SCHEME 1)
# We approximate the power-law kernel with exponential kernels with time constants τ
#   t = β .* 10 .^ (0:1/exp(1):3)
# and weights w where w = (1/τ)^ζ. ζ is an unknown adjustment factor that adjusts weights
# to better match the power-law kernel. 
function calc_optim_s1(
    β=1e-2; 
    dur=1e3β, 
    base=10.0, 
    start=0, 
    step=1/exp(1), 
    stop=3,
    init_mode="τ",
)
    # Determine time constants
    τ = β .* base .^ collect(start:step:stop)

    # Set initial parameters and bounds
    if init_mode == "τ"
        w = 1 ./ τ
    elseif init_mode == "τ + β"
        w = 1 ./ (τ .+ β)
    end

    # Synthesize time vector and pl kernel
    t = LinRange(log(1e-4), log(dur), 1000)
    kernel = log.(pl.(exp.(t), β))

    # Define function to compute approximation given τ and w
    approx(τ, w, ζ) = log.(β .* sum(pea_components(exp.(t), τ, w .^ ζ)))

    # Define loss function
    f = ζ -> loss(t, kernel, approx(τ, w, ζ[1]))

    # Compute optimization and return minimizer
    ζ̂ = Optim.minimizer(optimize(f, [0.8]; autodiff=:forward))[1]
    
    # Return results (for convenience, return τ, w, and ζ)
    τ, w .^ ζ̂, ζ̂
end

# calc_optim_s2 (SCHEME 2)
# We approximate the power-law kernel with exponential kernels with time constants τ
# and weights w. Time constants are forced to be:
#   τ = β .* 10 .^ (0:1/exp(1):3)
# but weights are free and numerically optimized.
function calc_optim_s2(
    β=1e-2; 
    dur=1e3β, 
    base=10.0, 
    start=0, 
    step=1/exp(1), 
    stop=3,
    init_mode="τ",
)
    # Determine time constants
    τ = β .* base .^ collect(start:step:stop)

    # Set initial parameters and bounds
    if init_mode == "τ"
        w = 1 ./ τ
    elseif init_mode == "τ + β"
        w = 1 ./ (τ .+ β)
    end

    # Synthesize time vector and pl kernel
    t = LinRange(log(1e-4), log(dur), 1000)
    kernel = log.(pl.(exp.(t), β))

    # Define function to compute approximation given τ and w
    approx(τ, w) = log.(β .* sum(pea_components(exp.(t), τ, exp.(w))))

    # Define loss function
    f = w -> loss(t, kernel, approx(τ, w))

    # Compute optimization and return minimizer
    ŵ = Optim.minimizer(optimize(f, log.(w); autodiff=:forward))
    
    # Return results (for convenience, return τ, w, and ζ)
    τ, exp.(ŵ)
end

# Calc heuristic approximation to the weights
function calc_heuristic_s1(
    β=1e-2; 
    base=10.0, 
    start=0, 
    step=1/exp(1), 
    stop=3,
)
    τ = β .* base .^ (start:step:stop)
    w = 1 ./ (τ .+ β)
    c = 2β * sum(w .* exp.(-β ./ τ))
    w = w ./ c
    return τ, w
end

# Figure 1
# Show target power-law kernel and approximation, as well as each component contributing
# to the fit, on a log-log scale.
function fig1(
    τ::Vector,
    w::Vector;
    β=1e-2, 
    fs=10e3, 
    xscale=log10, 
    yscale=log10,
    xlabel="Time (s)",
    ylabel="Amplitude (a.u.)",
    size=(1000, 700),
    dur=β*1e3,
    ylims=(1e-5, 5e0),
    xlims=(β*1e-1, β*1e3),
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
    plot_legend=false,
    plot_τ=false,
    plot_l=false,
    timevec_mode="log",
)
    # Create time vector & compute PEA
    if timevec_mode == "linear"
        t = timevec(dur, fs)
    else
        t = exp.(LinRange(log(1e-4), log(dur), 1000))
    end
    if length(τ) > 0
        g_comp = β .* pea_components(t, τ, w)
        g = β .* pea(t, τ, w)
    end
    f = pl.(t, β)

    # Select colors for each τᵢ
    if length(τ) == 1 || length(τ) == 0
        colors = [:black]
    else
        colors = get(colorschemes[colorscheme], LinRange(0.0, 1.0, length(τ)))
    end

    # Start figure with reference lines:
    #   1) horizontal gridline at 0.5
    # hlines!(ax, [0.5]; linestyle=:dash, color=:gray)
    if plot_τ vlines!(ax, τ; color=colors, linewidth=0.25) end

    # Plot results:
    #   1) Target in red
    #   2) Approximation components in colored lines, matching to legend
    #   3) Composite approximation in pink
    lines!(ax, t[2:end], f[2:end]; color=:red)
    if length(τ) > 0
        map(zip(g_comp, τ, colors)) do (gᵢ, τᵢ, colorᵢ)
            lines!(ax, t[2:end], gᵢ[2:end]; color=colorᵢ, label=string(round(τᵢ*1e3)) * " ms")
        end
    end
    if length(τ) > 1
        lines!(ax, t[2:end], g[2:end]; color=:gray, linestyle=:dash)
    end

    # Adjust limits
    ylims!(ax, ylims)
    xlims!(ax, xlims)

    # Add loss value to upper right
    l = loss(t, log.(f), log.(g))
    if plot_l
        text!(ax, [xlims[2]/1.1], [1.0]; text=string(round(l; digits=2)), align=(:right, :bottom), color=:gray)
    end

    # If size is small enough, skip every other xtick and ytick
    start = floor(log10(xlims[1]))
    stop = ceil(log10(xlims[2]))
    ticks = (start:1:stop)

#    if size[1] > 300
#        ax.xticks = 10 .^ ticks#, "10^−" .* string.(ticks))
#    end
    
    # Add legend
    if plot_legend axislegend(ax) end
    fig
end

# Convenient methods for fig1
fig1(β::Float64=1e-2; scheme=calc_optim_s1, kwargs...) = fig1(scheme(β)[1:2]...; β=β, kwargs...)
function fig1(β::Float64, ζ::Float64; base=10.0, start=0, step=1/exp(1), stop=3, kwargs...) 
    τ = β .* base .^ (start:step:stop)
    fig1(τ, 1 ./ τ .^ ζ; β=β, kwargs...)
end

# Figure 2
# Show Figure 1, but over multiple linear time scales
function fig2(args...; scheme=calc_optim_s1, β=1e-2, timescales=[(0.0, β), (0.0, 10*β), (0.0, 100*β)], kwargs...)
    # Fit parameters using calc_optim
    opt = scheme(β)

    # Create overall figure
    fig = Figure(; size=(1000, 300))

    # Loop and create each individual axis and plot over requested time scale
    for (idx, ts) in enumerate(timescales)
        ax = Axis(fig[1, idx])
        fig1(opt[1:2]...; β=β, dur=ts[2], fig=fig, ax=ax, xlims=ts, ylims=(0.0, 1.0))
    end
    fig
end

# Figure 3
# Show that optimization of ζ is well posed
function fig3(β=1e-2; base=10.0, start=0, step=1/exp(1), stop=3, fs=10e3, scale=log, yscale=identity)
    # Determine time constants
    τ = β .* base .^ collect(start:step:stop)

    # Set initial parameters and bounds
    w = 1 ./ τ

    # Synthesize time vector and pl kernel
    t = timevec(10.0, fs)
    kernel = pl.(t, β)

    # Define loss function (from scalar-valued x to scalar-valued loss)
    f = x -> loss(t, kernel, β .* sum(pea_components(t, τ, w .^ x)); scale=scale) 

    # Calculate losses
    Ζ = LogRange(2e-1, 2e0, 51)
    losses = map(f, Ζ)

    # Plot
    fig = Figure()
    ax = Axis(
        fig[1, 1]; 
        xscale=log10, 
        xminorticksvisible=true, 
        xminorticks=IntervalsBetween(9),
        yscale=yscale,
    )
    lines!(ax, Ζ, losses; color=:black)

    # Add vertical indicator for ζ estimate from `calc_optim`
    opt = calc_optim_s1(β; base=base, start=start, step=step, stop=stop, fs=fs, scale=scale)
    vlines!(ax, [opt[3]]; color=:red)
    
    # Add labels, adjust ticks, etc.
    ax.ylabel = "Loss (log-scaled)"
    ax.xlabel = "ζ (a.u.)"
    fig
end

# Figure 4
# Show relationship between β and ζ
function fig4(B=LogRange(1e-3, 1e-1, 15); scheme=calc_optim_s1, kwargs...)
    # Grab each estimate of ζ
    Z = map(x -> scheme(x; kwargs...)[3], B)

    # Create figure
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xscale=log10,
        yscale=log10,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        yminorticksvisible=true,
        yminorticks=IntervalsBetween(9),
    )

    # Plot
    scatter!(ax, B, Z)

    # Adjust labels, ticks, etc.
    ax.xlabel = "β (s)"
    ax.ylabel = "ζ (a.u.)"

    fig
end

# Figure 5
# Show relationship between β and loss metric
function fig5(B=LogRange(1e-3, 1e-1, 21); dur=1e1, fs=10e3, scheme=calc_optim_s1, kwargs...)
    # Grab loss for each β
    losses = map(B) do β 
        t = timevec(1000β, fs)
        k = pl.(t, β)
        loss(t, k, β .* pea(t, scheme(β; kwargs...)[1:2]...); scale=log10)
    end

    # Create figure
    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xscale=log10,
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
    )

    # Plot
    scatter!(ax, B, losses)

    # Adjust labels, ticks, etc.
    ax.xlabel = "β (s)"
    ax.ylabel = "ζ (a.u.)"

    fig
end

# Figure 6
# Show relationships for base, step, stop using numerical optimization
function fig6(; β=1e-2, dur=1e1)
    # Figure
    fig = Figure()

    # Effect of base
    steps = LogRange(1e-1, 1e0, 51)
    losses = pmap(steps) do step
        # Synthesize time vector and pl kernel
        t = LinRange(log(1e-4), log(dur), 1000)
        kernel = log.(pl.(exp.(t), β))

        # Define loss function
        τ, w = calc_optim_s4(β; base=10.0, step=step, stop=3.0)
        approx = log.(β .* sum(pea_components(exp.(t), τ, w)))
        loss(t, kernel, approx)
    end

    ax = Axis(
        fig[1, 1]; 
        xscale=log10, 
        yscale=log10, 
        yminorticksvisible=true, 
        yminorticks=IntervalsBetween(9)
    )
    lines!(ax, steps, losses; color=:black)
    ax.xlabel = "Exponent step size"
    ax.ylabel = "Log-log scale loss"
    step_min = steps[argmin(losses)]
    vlines!(ax, [step_min]; color=:red)
    text!(ax, [step_min*1.02], [1e1]; text="$(round(step_min; digits=3))", color=:red)
    fig
end

# Figure 7
# Show effects of maximum time for the loss function on the result of numerically optimized
# PLA weights
function fig7(β=1e-2; fs=10e3)
    # Determine settings
    maxdurs = LogRange(1e0, 1e4, 30)
    stops = 1:1:5

    # Compute time vector with log-spaced samples
    L = map(stops) do stop
        map(maxdurs) do maxdur
            t = LogRange(1/fs, maxdur, 1000)
            τ, w = calc_optim_s2(; stop=stop, dur=maxdur)
            k = log.(pl.(t, β))
            k̂ = log.(β .* sum(pea_components(t, τ, w)))
            loss(t, k, k̂)
        end
    end

    # Plot loss as figure
    fig = Figure()
    ax = Axis(fig[1, 1]; xscale=log10, yscale=log10)
    map(zip(L, stops)) do (ll, stop)
        lines!(ax, maxdurs, ll; label=string(stop))
    end

    axislegend()

    ylims!(ax, 1e-2, 1e7)
    xlims!(ax, 1e-1, 1e5)

    ax.xlabel = "Upper limit of time in loss function (s)"
    ax.ylabel = "Loss"

    fig
end

end # module PowerlawApproximation
