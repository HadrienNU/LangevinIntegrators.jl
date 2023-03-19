#=
timestep_harmonic.jl

Compute variance of position for an harmonic potential for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using Distributions

let
    force = ForceFromPotential("Harmonic")
    γ = 1.0
    β = 1.0
    t_range = LinRange(1e-3, 5e-1, 8)
    init_conds = Dict(
        "position" => Dict("type" => "gaussian", "std" => 1.0),
        "velocity" => Dict("type" => "gaussian", "std" => sqrt(1 / β)),
    )
    # inspectdr()
    plot(t_range, ones(length(t_range)), xaxis = :log, label = "Theory")
    for int_class in [ABOBA, BAOAB, BBK, GJF]

        var_position = zeros(size(t_range))
        err_var = zeros(size(t_range))
        for (n, Δt) in enumerate(t_range)
            params = TrajsParams(n_steps = 1e6, n_trajs = 250, n_save_iters = 5)
            println(String(Symbol(int_class)), " ", Δt)

            integrator = int_class(force, β, γ, 1.0, Δt, 1) # Change also initial conidition
            trajs =
                run_trajectories(integrator; params = params, init_conds_args = init_conds)
            x = vcat([trj.xt[1][:, 1] .^ 2 for trj in trajs]...)
            var_position[n] = mean(x)
            err_var[n] =
                quantile(TDist(length(x) - 1), 1 - 0.05 / 2) * std(x) / sqrt(length(x))
        end
        plot!(
            t_range,
            var_position,
            yerr = err_var,
            xaxis = :log,
            marker = (:circle, 1),
            label = String(Symbol(int_class)),
        )
    end
    xlabel!("Δt")
    ylabel!("<x_t^2>")
end
