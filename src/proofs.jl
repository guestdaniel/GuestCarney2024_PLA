# proofs.jl
#
# Various simple figures showing graphical equivalents of certain useful proofs

export proof_decay_and_cutoff, proof_decay_vs_tau, proof_approximate_vs_exact_ω, proof_filt_vs_direct, proof_match_kernel_match_adapt

# Useful identities/relations
τ_to_f(τ) = 1/(2π * τ)
f_to_τ(f) = 1/(2π * f)
f_to_ω(f, fs) = 2π * f / fs
ω_to_d(ω) = cos(ω) - 1 + sqrt(cos(ω)^2 - 4*cos(ω) + 3)
τ_to_d(τ, fs) = ω_to_d(f_to_ω(τ_to_f(τ), fs))
τ_to_d_quick(τ, fs) = 1 - exp(-1/fs / τ)

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
    xlims!(ax, 0.0, 5τ)
    ax.title = string(d)
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

# proof_filt_vs_direct()
function proof_filt_vs_direct(τ=1e-2; fs=10e3)
    # Calc time vector 
    t = timevec(1.0, fs)

    # Convert τ to decay coefficients
    d = τ_to_d_quick(τ, fs)
    b = [1]
    a = [1, -(1-d)]
    x = vcat(1.0, zeros(Int(fs)-1))
    y1 = filt(b, a, x)

    # Do direct
    y2 = zeros(size(y1))
    for i in eachindex(y1)
        if i == 1
            y2[i] = x[i]
        else
            y2[i] = x[i] + (1-d)*y2[i-1]
        end
    end

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, t, y1)
    lines!(ax, t, y2; linestyle=:dash)
    xlims!(ax, 0.0, 5τ)
    vlines!(ax, [τ])
    hlines!(ax, [1/exp(1)])
    fig
end

# proof_match_kernel_match_adapt
function proof_match_kernel_match_adapt(β=5e-4; fs=10e3)
    # Create time axis 
    t = timevec(1.0, fs)

    # Select automatically perfect inputs
    τ, w = calc_optim_s2(β)
    d = τ_to_d_quick.(τ, fs)

    # Input: constant unit amplitude
    x = ones(samples(1.0, fs))

    # # Algorithm #1: PLA
    I = zeros(length(x))
    y_pla = zeros(length(x))
    for i in eachindex(y_pla)
        I[i] = sum(y_pla[1:i] ./ (t[i] .- t[1:i] .+ β)) * (1/fs)
        y_pla[i] = max(x[i] - I[i], 0.0)
    end
    I_pla = I
    I_pla_conv = conv(y_pla, 1 ./ (t .+ β))

    # Algorithm #2: PEA
    I = 0.0
    I_pea = zeros(size(x))
    E = zeros(length(τ))
    y_pea = zeros(length(x))
    for i in eachindex(y_pea)
        y_pea[i] = max(x[i] - I, 0.0)
        for j = 1:length(τ)
            if i == 1
                E[j] = y_pea[i]
            else
                E[j] = y_pea[i] + (1-d[j]) * E[j]
            end
        end
        I = 1/fs * sum(w .* E)
        I_pea[i] = I
    end

    # Plot
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, t, I_pla)
    lines!(ax, t, I_pla_conv[1:Int(fs)] ./ fs; color=:red, linestyle=:dash)
    lines!(ax, t, I_pea)
    xlims!(ax, 0.0, 5β)
    ax = Axis(fig[1, 2])
    lines!(ax, t, y_pla)
    lines!(ax, t, y_pea)
    xlims!(ax, 0.0, 5β)
    fig
end