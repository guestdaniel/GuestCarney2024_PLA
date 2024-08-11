module PowerlawApproximation

using AuditorySignalUtils
using CairoMakie
using Optim
using Colors
using ColorSchemes

export pl, e, pea, fig1, fig2, loss, calc_optim

# Compute power-law function `pl` or exponential function `e`
pl(t, β) = β / (t + β)
e(t, τ) = exp(-t/τ)

# Compute exponential to approximate power-law 
pea(t, τ) = map(τᵢ -> 1/τᵢ * e.(t, τᵢ), τ)

# Compute simple squared error loss
function loss(x, y; fs=10e3)
    sum((y .- x) .^ 2 .* (1/fs))
end

# Optimize single time constant
function calc_optim(β; base=10.0, start=0, step=1, stop=2, fs=10e3)
    # Determine time constants
    τ = β .* base .^ collect(start:step:stop)

    # Set initial parameters and bounds
    ζ = 0.2

    # Synthesize time vector and pl kernel
    t = timevec(10.0, fs)
    kernel = pl.(t, β)

    # Define loss function
    f = x -> loss(kernel, β * sum(τ .^ x[1] .* pea(t, τ)))  # TODO: try log-scaled 

    # Compute minimization
    τ, Optim.minimizer(optimize(f, [ζ]))
end

# Plot PEA fit to kernel
function fig1(; 
    β=1e-2, 
    base=10.0,
    start=0,
    step=1/exp(1),
    stop=3,
    kwargs...
)
    τ = β .* base .^ (start:step:stop)
    fig1(τ; β=β, kwargs...)
end

function fig1(
    τ,
    ζ=1.0,
    c=τ .^ ζ;
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
    plot_kernel_pl=true,
    plot_kernel_pea=true,
    plot_kernel_pea_components=true,
)
    # Create time vector & compute PEA
    t = timevec(dur, fs)
    g = β .* pea(t, τ)

    # Select colors for each τᵢ
    colors = get(colorschemes[colorscheme], LinRange(0.0, 1.0, length(τ)))

    # Start figure with reference lines:
    #   1) horizontal gridline at 0.5
    hlines!(ax, [0.5]; linestyle=:dash, color=:gray)

    # Optionally, plot each component of the PEA
    if plot_kernel_pea_components   
        map(zip(c .* g, τ, colors)) do (gᵢ, τᵢ, colorᵢ)
            lines!(ax, t[2:end], gᵢ[2:end]; color=colorᵢ, label=string(round(τᵢ*1e3)) * " ms")
        end
    end

    # Optionally, plot composite PEA kernel
    if plot_kernel_pea
        lines!(ax, t[2:end], sum(c .* g)[2:end]; color=:pink, linestyle=:dash)
    end

    # Optionally, plot PL kernel
    if plot_kernel_pl
        lines!(ax, t[2:end], pl.(t[2:end], β); color=:red)
    end

    # Adjust limits
    ylims!(ax, ylims)
    xlims!(ax, xlims)
    
    # Add legend
    if plot_kernel_pea_components axislegend(ax) end
    fig
end

# Plot PEA fit to kernel
function fig1(; 
    β=1e-2, 
    base=10.0,
    start=0,
    step=1/exp(1),
    stop=3,
    kwargs...
)
    τ = β .* base .^ (start:step:stop)
    fig1(τ; β=β, kwargs...)
end

function fig2(
    τ,
    ζ=1.0,
    c=τ .^ ζ;
    β=1e-2, 
    fs=10e3, 
    xscale=log10, 
    yscale=log10,
    xlabel="Time (s)",
    ylabel="Amplitude (a.u.)",
    size=(300, 1500),
    dur=1e1,
    ylims=(1e-5, 5e0),
    xlims=(1/fs / 10, dur * 10),
    fig=Figure(; size=size), 
    colorscheme=:glasgow,
    plot_kernel_pl=true,
    plot_kernel_pea=true,
    plot_kernel_pea_components=true,
)
    # Create time vector & compute PEA
    t = timevec(dur, fs)
    g = β .* pea(t, τ)

    # Select colors for each τᵢ
    colors = get(colorschemes[colorscheme], LinRange(0.0, 1.0, length(τ)))

    for i in eachindex(τ)
        # Create axis
        ax=Axis(
            fig[i, 1]; 
            xscale=xscale,
            yscale=yscale,
            xlabel=xlabel,
            ylabel=ylabel,
            xminorticksvisible=true,
            xminorticks=IntervalsBetween(9),
            yminorticksvisible=true,
            yminorticks=IntervalsBetween(9),
        )
        # Start figure with reference lines:
        #   1) horizontal gridline at 0.5
        #   2) vertical gridline at β
        #   3) vertical gridline at τᵢ
        hlines!(ax, [0.5]; color=:gray)
        vlines!(ax, [β]; color=:red)
        vlines!(ax, [τ[i]]; color=colors[i])

        # Optionally, plot each component of the PEA
        if plot_kernel_pea_components
            lines!(ax, t[2:end], c[i] .* g[i][2:end]; color=colors[i], label=string(round(τ[i]*1e3)) * " ms")
        end

        # Optionally, plot composite PEA kernel
        if plot_kernel_pea
            lines!(ax, t[2:end], sum(c[1:i] .* g[1:i])[2:end]; color=:pink, linestyle=:dash)
        end

        # Optionally, plot PL kernel
        if plot_kernel_pl
            lines!(ax, t[2:end], pl.(t[2:end], β); color=:red)
        end

        # Hide extraneous stuff
        if i < length(τ)
            hidexdecorations!(ax, ticks=false, grid=false)
        end

        # Adjust limits
        ylims!(ax, ylims)
        xlims!(ax, xlims)
    end
    fig
end

end # module PowerlawApproximation
