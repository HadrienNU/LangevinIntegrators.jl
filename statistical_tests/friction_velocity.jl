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
function velocity_diff(xt, vt, lambda, Δt)
    ut = lambda*(xt[3:end]-xt[2:end-1])/Δt+ (1.0 -lambda)*(xt[2:end-1]-xt[1:end-2])/Δt
    return xt[3:end-1], vt[2:end-2]-ut[1:end-1]
    # return vt[3:end-1]- ut[2:end], vt[2:end-2]-ut[1:end-1]
end


let
    force=ForceFromPotential("Flat")
    γ_range = LinRange(0.1,10.0,10)
    β = 1.0
    t_range = LinRange(1e-5,5e-2,5)
    lambda_range = LinRange(-0.5,1.0,7)
    init_conds=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"gaussian", "std"=> sqrt(1/β)))
    inspectdr()

    plot(γ_range, γ_range./3, label = "Theory")
    plot!(γ_range, -γ_range./12, label = "Theory")
    int_class = VEC

    Δt=1e-3
    for int_class in [BAOAB,OBABO,BBK,GJF,VEC] #
        var_meanvel=zeros(length(γ_range),2)
        err_meanvel=zeros(length(γ_range),2)
        var_stdvel=zeros(length(γ_range),2)
        for (n,γ) in enumerate(γ_range)
            params=TrajsParams(n_steps = 1e5, n_trajs = 15, n_save_iters = 1)
            println(String(Symbol(int_class))," ", Δt)
            integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
            trajs=run_trajectories(integrator; params = params,to_save=["x","v","v_mid"], init_conds_args=init_conds)
            x = vcat([hcat(velocity_diff(trj.xt[1][1:end,1],trj.xt[2][1:end,1],0.5, Δt)...) for trj in trajs]...)
            var_meanvel[n,:] = mean(x,dims=1)
            cov_mat = cov(x)
            var_stdvel[n,1]=cov_mat[1,1]/(Δt)
            var_stdvel[n,2]=cov_mat[1,2]/(Δt)
            err_meanvel[n,:] =quantile(TDist(size(x)[1]-1), 1 - 0.05/2) * std(x,dims=1)/sqrt(size(x)[1])
        end
         # plot!(lambda_range, var_meanvel[:,1], marker=(:circle,10),label="$(String(Symbol(int_class))) $Δt" )
         plot!(γ_range, var_stdvel[:,2], marker=(:auto,10),label="$(String(Symbol(int_class))) std $Δt" )
    end
    xlabel!("γ")
    ylabel!("<vₜ-u^λₜ vₜ-u^λₜ>")
end
