#=
precompile_all.jl

Run all integrators for a small number of time step for
=#

using LangevinIntegrators

let
    force=ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    n_step = 10
    Δt = 1e-2

    for int_class in [ABOBA,BAOAB,OBABO,BBK,GJF,VEC]
        params=TrajsParams(n_steps = n_step, n_trajs = 2, n_save_iters = n_step)
        println(String(Symbol(int_class)))
        integrator=int_class(force, β , γ, 1.0, Δt, 1)
        trajs=run_trajectories(integrator; params = params)
    end
end
