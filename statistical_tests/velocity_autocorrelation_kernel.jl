#=
main.jl

Compute velocity autocorrelation function for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using StatsBase

let
    force = ForceFromPotential("Flat")

    params = TrajsParams(n_steps = 5 * 10^5, n_trajs = 50, n_save_iters = 1)
    Δt = 1e-3
    factor_length = 50
    tps = LinRange(
        0,
        div(params.n_steps, factor_length) * Δt,
        div(params.n_steps, factor_length) + 1,
    )
    α = 5.0
    max_τ = 2.0
    max_step = trunc(Int, max_τ / Δt)
    kernel = 20.0 * exp.(-α * LinRange(0, max_step * Δt, max_step + 1))
    β = 1.0
    Δ = sqrt(4 * kernel[1] - α^2)
    init_conds = Dict(
        "position" => Dict("type" => "gaussian", "std" => 1.0),
        "velocity" => Dict("type" => "gaussian", "std" => sqrt(1 / β)),
    )
    inspectdr()
    plot(
        tps,
        exp.(-0.5 * α * tps) .* (cos.(0.5 * Δ * tps) + α / Δ * sin.(0.5 * Δ * tps)) / β,
        label = "Theory",
    )
    for int_class in [BBK_Kernel, GJF_Kernel]
        println(String(Symbol(int_class)))
        integrator = int_class(force, β, kernel, 1.0, Δt, 1)

        trajs = run_trajectories(
            integrator;
            params = params,
            to_save = ["v"],
            init_conds_args = init_conds,
        )
        #Compute velocity auto correlation and plot it
        corr = zeros(div(params.n_steps, factor_length) + 1)
        for trj in trajs
            corr .+=
                correlation(
                    trj.xt[1][10^3:end, 0],
                    trunc = div(params.n_steps, factor_length) + 1,
                ) / params.n_trajs
            # corr.+=autocov(trj.xt[1][10^3:end],0:div(params.n_steps,factor_length))[:,1]/params.n_trajs
        end

        plot!(tps, corr, label = String(Symbol(int_class)))
    end
    xlabel!("t")
    ylabel!("<v_0,v_t>")

end
