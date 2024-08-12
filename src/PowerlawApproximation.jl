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

# Compute exponential to approximate power-law 
pea(t, τ, w=ones(length(τ))) = sum(map(x -> x[2] .* e.(t, x[1]), zip(τ, w)))
pea_components(t, τ, w=ones(length(τ))) = map(x -> x[2] .* e.(t, x[1]), zip(τ, w))

# Compute simple squared error loss
function loss(t, x, y)
    trapz((t), (y .- x) .^ 2)
end

# Optimize single time constant
function calc_optim(β; base=10.0, start=0, step=1/exp(1), stop=3, fs=10e3)
    # Determine time constants
    τ = β .* base .^ collect(start:step:stop)

    # Set initial parameters and bounds
    ζ = 0.2

    # Synthesize time vector and pl kernel
    t = timevec(10.0, fs)
    kernel = pl.(t, β)

    # Define loss function
    f = x -> loss(t, kernel, β .* sum((1 ./ τ.^2) .^ x[1] .* pea_components(t, τ))) 

    # Compute minimization
    ζ̂ = Optim.minimizer(optimize(f, [ζ]))
    τ, (1 ./ τ .^ 2) .^ ζ̂, ζ̂
end

function fig1(
    τ,
    w,
    ζ;
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
        #ylabel=ylabel,
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
    g_comp = β .* w .* pea_components(t, τ)
    g = β .* pea(t, τ, w)
    f = pl.(t, β)

    # Select colors for each τᵢ
    colors = get(colorschemes[colorscheme], LinRange(0.0, 1.0, length(τ)))

    # Start figure with reference lines:
    #   1) horizontal gridline at 0.5
    hlines!(ax, [0.5]; linestyle=:dash, color=:gray)

    # Optionally, plot each component of the PEA
    if plot_kernel_pea_components   
        map(zip(g_comp, τ, colors)) do (gᵢ, τᵢ, colorᵢ)
            lines!(ax, t[2:end], gᵢ[2:end]; color=colorᵢ, label=string(round(τᵢ*1e3)) * " ms")
        end
    end

    # Optionally, plot composite PEA kernel
    if plot_kernel_pea
        lines!(ax, t[2:end], g[2:end]; color=:pink, linestyle=:dash)
    end

    # Optionally, plot PL kernel
    if plot_kernel_pl
        lines!(ax, t[2:end], f[2:end]; color=:red)
    end

    # Adjust limits
    ylims!(ax, ylims)
    xlims!(ax, xlims)
    
    # Add legend
    if plot_kernel_pea_components axislegend(ax) end
    fig
end

# function fig2(
#     τ,
#     ζ=1.0,
#     c=τ .^ ζ;
#     β=1e-2, 
#     fs=10e3, 
#     xscale=log10, 
#     yscale=log10,
#     xlabel="Time (s)",
#     ylabel="Amplitude (a.u.)",
#     size=(360, 1000),
#     dur=1e1,
#     ylims=(yscale == identity ? 0.01 : 1e-3, yscale == identity ? 1.1e0 : 5e0),
#     xlims=(xscale == identity ? 0.0 : 1/fs / 10, xscale == identity ? β*3e1 : dur * 10),
#     fig=Figure(; size=size), 
#     colorscheme=:glasgow,
#     base=10.0,
# )
#     # Create time vector & compute PEA
#     t = timevec(dur, fs)
#     g = β .* pea(t, τ)

#     # Select colors for each τᵢ
#     colors = get(colorschemes[colorscheme], LinRange(0.0, 1.0, length(τ)))

#     for i in eachindex(τ)
#         # Create axis
#         ax=Axis(
#             fig[i, 1]; 
#             xscale=xscale,
#             yscale=yscale,
#             xlabel=xlabel,
#             #ylabel=ylabel,
#             xminorticksvisible=true,
#             xminorticks=IntervalsBetween(9),
#             yminorticksvisible=true,
#             yminorticks=IntervalsBetween(9),
#         )
#         # Start figure with reference lines:
#         #   1) horizontal gridline at 0.5
#         #   2) vertical gridline at β
#         #   3) vertical gridline at τᵢ
#         # hlines!(ax, [0.5]; color=:gray)
#         vlines!(ax, [β]; color=:black, linewidth=0.5)
#         vlines!(ax, [τ[i]]; color=colors[i])

#         # Compute each component of the PEA
#         kernel_pl = pl.(t[2:end], β)
#         gᵢ = c[i] .* g[i][2:end]
#         ĝ = sum(c[1:i] .* g[1:i])[2:end]

#         # Optionally, plot each component of the PEA
#         t_highlight = (t .> (τ[i] ^ base)) .& (t .> (τ[i] ^ -base))
#         lines!(ax, t[2:end], kernel_pl; color=:black, linewidth=0.5)
#         lines!(ax, t[2:end], gᵢ; color=colors[i], label=string(round(τ[i]*1e3)) * " ms")
#         lines!(ax, t[2:end], ĝ; color=:pink)

#         # Hide extraneous stuff
#         if i < length(τ)
#             hidexdecorations!(ax, ticks=false, grid=false, minorticks=false)
#         end

#         # Adjust limits
#         ylims!(ax, ylims)
#         xlims!(ax, xlims)

#         # Text
#         text!(ax, [1e2], [1.0]; text="$(round(1e5*loss(kernel_pl, ĝ; fs=fs); digits=5))", align=(:right, :bottom))

#         # Adjust spacing
#         rowgap!(fig.layout, Relative(2e-3))
#     end
#     fig
# end

end # module PowerlawApproximation
