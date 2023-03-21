#=
main.jl

Compute diffusion coefficients for various set of integrators and compare it to expected value
=#

using Plots
using StatsBase
using LangevinIntegrators

let
    force = ForceFromPotential("Harmonic")
    γ = 1.0
    β = 1.0
    Δt = 5e-3
    x_range = LinRange(-4.0, 4.0, 100)
    plot(x_range, 0.5 * (x_range) .^ 2, label = "Theory")
    for int_class in [ABOBA, BAOAB, BBK, GJF]
        params = TrajsParams(n_steps = 1000000, n_trajs = 1, n_save_iters = 5)
        println(String(Symbol(int_class)))
        integrator = int_class(force, β, γ, 1.0, Δt, 1) # Change also initial conidition
        trajs = run_trajectories(integrator; params = params)
        h = fit(Histogram, trajs[1].xt[1][:, 1], nbins = 50)
        log_hist = -log.(h.weights)
        plot!(
            0.5 * (h.edges[1][1:end-1] + h.edges[1][2:end]),
            log_hist .- minimum(log_hist),
            marker = (:circle, 1),
            label = String(Symbol(int_class)),
        )
        # plot!(xrange, histogram,label=String(Symbol(int_class)))
    end

    xlabel!("x")
    ylabel!("-log P(x)")
end
