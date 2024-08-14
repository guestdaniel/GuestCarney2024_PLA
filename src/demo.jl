# demo.jl
# 
# Nice miniature examples to show during lab meeting / grad meeting to demo the ideas

# Set theme
set_theme!()

# Example 0: Lab meeting
β = 1e-2
fig1([], []; xscale=identity, yscale=identity, xlims=(0.0, 5β), ylims=(0.0, 1.0))
fig1([β/5], [1/β]; xscale=identity, yscale=identity, xlims=(0.0, 5β), ylims=(0.0, 1.0))
fig1([β*5], [1/β]; xscale=identity, yscale=identity, xlims=(0.0, 5β), ylims=(0.0, 1.0))
fig1([β], [1/β]; xscale=identity, yscale=identity, xlims=(0.0, 5β), ylims=(0.0, 1.0))
fig1([β], [1/β]; xscale=identity, yscale=identity, xlims=(0.0, 500β), ylims=(0.0, 1.0))
fig1([β], [1/β])
fig1([10*β], [1/β])

τ = β .* 10 .^ (0.0:0.5:3.0)
fig1(τ, 10.0 .* ones(length(τ)))
fig1(τ, 10.0 .* ones(length(τ)); xscale=identity, yscale=identity, xlims=(0.0, 5β), ylims=(0.0, 2.0))
fig1(τ, 10.0 .* ones(length(τ)))
fig1(τ, 1.0 ./ τ)
fig1(τ, (1 ./ τ) .^ 0.9; xscale=identity, yscale=identity, xlims=(0.0, 500β), ylims=(0.0, 2.0))

# Example 1:
# Demonstrate basic philosophy of the modeling approach as a function of β
# Try: β=1e-3, β=1e-2, β=1e-1
β = 1e-2
τ = β .* 10.0 .^ (0.0:1/exp(1):3.0)
fig1(τ, 1 ./ τ .* LinRange(0.5, 1.0, length(τ)); β=β)
fig1(τ, 1 ./ τ .* LinRange(0.5, 1.0, length(τ)); β=β, xscale=identity, yscale=identity, xlims=(0.0, 10β), ylims=(0.0, 1.0))

# Example 2:
# Demonstrate philosophy of different ways of adjusting w
# Try: β=1e-3, β=1e-2, β=1e-1
β = 1e-2
τ = β .* 10.0 .^ (0.0:1/exp(1):3.0)
fig1(τ, 1 ./ τ .* LinRange(0.5, 1.0, length(τ)); β=β)  # linear ramp from 0.5 to 1
fig1(τ, 1 ./ τ .^ 0.83; β=β)                           # raised to power 

# Example 3: More complex strategy for adjusting w
fig1(τ, 1 ./ (τ .+ β))                             # use 1/(τ+β)

τ = β .* 10.0 .^ (0.0:1/exp(1):4.0)
w = 1 ./ τ
tar_at_0 = 1.0
val_at_0 = β * sum(w .* exp.(0 ./ τ))
r = tar_at_0 / val_at_0
fig1(τ, r .* w .* LinRange(0.5, 1.0, length(w)))