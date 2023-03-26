#=
timestep_harmonic_velocity.jl

Compute variance of velocity for an harmonic potential for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators
using Distributions

"""
Compute on one traj the wanted observable
"""

function likelihood_vec_finite(xt, lambda, Δt, γ)
    ut = lambda * (xt[3:end] - xt[2:end-1]) / Δt + (1.0 - lambda) * (xt[2:end-1] - xt[1:end-2]) / Δt
    return likelihood_vec(xt[2:end-1], ut, Δt, γ)
end

function likelihood_vec(xt, vt, Δt, γ)
    Mqq = 2*(Δt^3)*γ/3
    Mqv = γ*(Δt^2)*(1-γ*Δt)
    Mvv = 2*γ*Δt*(1-γ*Δt+2*(Δt^2)*γ^2/3)
    detM = Mqq*Mvv-Mqv^2
    Mq = xt .+ vt.*Δt*(1-0.5*γ*Δt)
    Mv = vt.*(1-γ*Δt+0.5*(γ*Δt)^2)
    return Mqq/detM .*(vt[2:end].-Mv[1:end-1]).^2 #Mvv/detM .*(xt[2:end].-Mq[1:end-1]).^2 #
end

let
    force = ForceFromPotential("Flat")
    γ_range = LinRange(0.01, 1.5, 10)
    β = 1.0
    t_range = LinRange(5e-4, 5e-2, 5)
    lambda_range = LinRange(-0.5, 1.0, 7)
    init_conds = Dict(
        "position" => Dict("type" => "Cste"),
        "velocity" => Dict("type" => "gaussian", "std" => sqrt(1 / β)),
    )
    inspectdr()
    plot()

    every = 1

    plot()
    # int_class = BBK
    # for (m,Δt) in enumerate(t_range)
    Δt = 1e-3
    for int_class in [VEC] #BAOAB, OBABO, BBK, GJF,
        var_finite = zeros(size(γ_range))
        var_vel = zeros(size(γ_range))
        for (n,γ) in enumerate(γ_range)
            params = TrajsParams(n_steps = 1e5, n_trajs = 15, n_save_iters = 1)
            println(String(Symbol(int_class)), " ", γ)
            integrator = int_class(force, β, γ, 1.0, Δt, 1) # Change also initial conidition
            trajs = run_trajectories(
                integrator;
                params = params,
                to_save = ["x", "v", "v_mid"],
                init_conds_args = init_conds,
            )

            c2 = (1.0 - integrator.c₂) / Δt
            x_finite = vcat(
                [
                    likelihood_vec_finite(
                        trj.xt[1][:, 1],
                        0.5,
                        every * Δt,
                        c2,
                    ) for trj in trajs
                ]...,
            )
            println(x_finite[1:5])
            var_finite[n] = mean(x_finite)
            x = vcat(
                [
                    likelihood_vec(
                        trj.xt[1][:, 1],
                        trj.xt[2][:, 1],
                        every * Δt,
                        c2,
                    ) for trj in trajs
                ]...,
            )
            println(x[1:5])
            var_vel[n] = mean(x)

        end
        plot!(
            γ_range,
            var_vel,
            marker = (:auto, 10),
            label = "$(String(Symbol(int_class))) $Δt",
        )
        plot!(
            γ_range,
            var_finite,
            marker = (:auto, 10),
            label = "Finite $(String(Symbol(int_class))) $Δt",
        )
    end
    xlabel!("λ")
    ylabel!("L")
end
