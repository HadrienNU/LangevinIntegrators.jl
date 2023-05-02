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


let
    force=ForceFromPotential("Flat")
    γ_range = LinRange(0.1,10.0,10)
    β = 1.0
    t_range = LinRange(1e-5,5e-2,5)
    lambda_range = LinRange(-0.5,1.0,4)
    init_conds=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"gaussian", "std"=> sqrt(1/β)))
    inspectdr()

    plot()
    int_class = VEC

    # Δt=1e-2
    γ=1.0
    for int_class in [VEC] #
        var_meanvel=zeros(length(t_range),3)
        err_meanvel=zeros(length(t_range),3)
        # for (n,γ) in enumerate(γ_range)
        for (n,Δt) in enumerate(t_range)
            params=TrajsParams(n_steps = 1e5, n_trajs = 1000, n_save_iters = 1)
            println(String(Symbol(int_class))," ", Δt)
            integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
            trajs=run_trajectories(integrator; params = params,to_save=["x","ξ","ξ₂"], init_conds_args=init_conds)
            x = vcat([trj.xt[1][2:end,1].*trj.xt[2][1:end-1,1] for trj in trajs]...)
            var_meanvel[n,1] = mean(x)
            err_meanvel[n,1] =quantile(TDist(size(x)[1]-1), 1 - 0.05/2) * std(x)/sqrt(size(x)[1])
            x = vcat([trj.xt[1][1:end,1].*trj.xt[2][1:end,1] for trj in trajs]...)
            var_meanvel[n,2] = mean(x)
            err_meanvel[n,2] =quantile(TDist(size(x)[1]-1), 1 - 0.05/2) * std(x)/sqrt(size(x)[1])

            x = vcat([trj.xt[1][1:end-1,1].*trj.xt[2][2:end,1] for trj in trajs]...)
            var_meanvel[n,3] = mean(x)
            err_meanvel[n,3] =quantile(TDist(size(x)[1]-1), 1 - 0.05/2) * std(x)/sqrt(size(x)[1])


        end
         # plot!(lambda_range, var_meanvel[:,1], marker=(:circle,10),label="$(String(Symbol(int_class))) $Δt" )
         plot!(t_range, var_meanvel, yerr = err_meanvel, marker=(:auto,10),label="$(String(Symbol(int_class)))" )
    end
    xlabel!("γ")
    ylabel!("<vₜ-u^λₜ vₜ-u^λₜ>")
end
