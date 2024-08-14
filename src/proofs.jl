# proofs.jl
#
# Various simple figures showing graphical equivalents of certain useful proofs

export proof_decay_and_cutoff, proof_decay_vs_tau, proof_approximate_vs_exact_ω

# Useful identities/relations
τ_to_f(τ) = 1/(2π * τ)
f_to_τ(f) = 1/(2π * f)
f_to_ω(f, fs) = 2π * f / fs
ω_to_d(ω) = cos(ω) - 1 + sqrt(cos(ω)^2 - 4*cos(ω) + 3)
τ_to_d(τ, fs) = ω_to_d(f_to_ω(τ_to_f(τ), fs))

# proof_decay_and_cutoff()
# Visual proof of relationship between decay coefficient and cutoff frequency
function proof_decay_and_cutoff(ω=0.5*π; fs=10e3)
    # A single-pole IIR filter is defined by the difference equation:
    #   y[n] = d*x[n] + (1-d)*y[n-1]
    # Thus, it has the transfer function in the z-domain of
    #	H(z) = d / (1 - (1-d) * z^-1)
    # It can be shown that an exact solution relating the decay coefficient d
    # to the cutoff frequency (omega_c) is:
    #	d = -y + sqrt(y^2 + 2y) where y = 1 - cos(omega_c)

    # Calculate f and d from ω
    # f = ω*fs/(2*pi)
    d = ω_to_d(ω)
    b = [d]
    a = [1, -(1-d)]

    # Filter impulse calculate spectrum of impulse repsonse
    unit_impulse = vcat(1.0, zeros(Int(fs)-1))
    ir = filt(b, a, unit_impulse)

    # Create figure
    fig = Figure()
    ax = Axis(fig[1, 1])
    Ω = LinRange(0, 2π, length(ir))
    H = 20.0 .* log10.(abs.(fft(ir)))
    lines!(ax, Ω, H)
    ax.xlabel = "Normalized frequency (x π rad/sample)"
    ax.ylabel = "Magnitude (dB)"
    vlines!(ax, [ω])
    hlines!(ax, [-3])
    fig
end

# proof_decay_vs_tau()
# Prove relationship between decay vs tau
function proof_decay_vs_tau(τ; fs=10e3)
    d = τ_to_d(τ, fs)
    b = [1]
    a = [1, -(1-d)]
    y = filt(b, a, vcat(1.0, zeros(Int(fs)-1)))

    fig = Figure()
    ax = Axis(fig[1, 1])
    t = timevec(1.0, fs)
    lines!(ax, t, y)
    vlines!(ax, [τ])
    hlines!(ax, [1/exp(1)])
    fig
end

# proof_approximate_vs_exact_ω()
function proof_approximate_vs_exact_ω(; fs=100e3)
    T = LogRange(1e2, 1e3, 101)
    Ω = f_to_ω.(τ_to_f.(T), fs)
    D_exact = ω_to_d.(Ω)
    D_approximate = exp.(-(1/fs) ./ T)

    fig = Figure()
    ax = Axis(fig[1, 1])#, xscale=log10, yscale=log10)
    lines!(ax, T, D_exact; color=:blue)
    lines!(ax, T, 1 .- D_approximate; color=:gray, linestyle=:dash)
    ax.xlabel = "τ (s)"
    ax.ylabel = "d (a.u.)"
    fig
end