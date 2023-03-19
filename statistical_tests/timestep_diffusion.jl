#=
main.jl

Compute diffusion coefficients for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using Distributions

let
    force = ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    n_step = 1e5
    t_range = LinRange(5e-4, 0.2, 5)
    init_conds = Dict(
        "position" => Dict("type" => "Cste"),
        "velocity" => Dict("type" => "gaussian", "std" => sqrt(1 / β)),
    )
    inspectdr()
    plot(t_range, ones(length(t_range)), xaxis = :log, label = "Theory")
    # for int_class in [ABOBA,BAOAB,BBK,GJF]
    #
    #     diff_coeff=zeros(size(t_range))
    #     err_diff = zeros(size(t_range))
    #     for (n,Δt) in enumerate(t_range)
    #         params=TrajsParams(n_steps = n_step, n_trajs = 1e4, n_save_iters = n_step)
    #         println(String(Symbol(int_class))," ", Δt)
    #
    #         integrator=int_class(force, β , γ, 1.0, Δt, 1)
    #         trajs=run_trajectories(integrator; params = params, init_conds_args=init_conds)
    #         end_points = [(trj.xt[1][end,1])^2/(2*params.n_steps*Δt) for trj in trajs]
    #         diff_coeff[n] = mean(end_points)
    #         err_diff[n] =quantile(TDist(params.n_trajs-1), 1 - 0.05/2) * std(end_points)/sqrt(params.n_trajs)
    #     end
    #      plot!(t_range, diff_coeff, yerr = err_diff,marker=(:circle,1),label=String(Symbol(int_class)))  # L'erreur est alors sqrt(var_diff)/params.n_trajs
    # end

    # Kernel integrator
    α = 20.0
    for int_class in [BBK_Kernel, GJF_Kernel]
        diff_coeff = zeros(size(t_range))
        err_diff = zeros(size(t_range))
        for (n, Δt) in enumerate(t_range)
            params = TrajsParams(n_steps = n_step, n_trajs = 1e4, n_save_iters = n_step)
            println(String(Symbol(int_class)), " ", Δt)

            max_step = trunc(Int, 2.0 / Δt)
            γ = 20.0 * exp.(-α * LinRange(0, max_step * Δt, max_step + 1))

            integrator = int_class(force, β, γ, 1.0, Δt, 1)
            trajs =
                run_trajectories(integrator; params = params, init_conds_args = init_conds)
            end_points =
                [(trj.xt[1][end, 1])^2 / (2 * params.n_steps * Δt) for trj in trajs]
            diff_coeff[n] = mean(end_points)
            err_diff[n] =
                quantile(TDist(params.n_trajs - 1), 1 - 0.05 / 2) * std(end_points) /
                sqrt(params.n_trajs)
        end
        plot!(
            t_range,
            diff_coeff,
            yerr = err_diff,
            marker = (:circle, 1),
            label = String(Symbol(int_class)),
        )  # L'erreur est alors sqrt(var_diff)/params.n_trajs
    end

    xlabel!("Δt")
    ylabel!("<(x_t-x_0)^2>")
end
