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
    return vt[3:end-1]- ut[2:end], vt[2:end-2]-ut[1:end-1]
end


let
    force=ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    t_range = LinRange(5e-4,5e-2,5)
    lambda_range = LinRange(-0.5,1.0,7)
    init_conds=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"gaussian", "std"=> sqrt(1/β)))
    inspectdr()

    plot(lambda_range, 2*(lambda_range.^2 + (1.0 .-lambda_range).^2)./3, label = "Theory")
    plot!(lambda_range, -2 .*lambda_range.*(1.0 .-lambda_range)./6, label = "Theory")
    # int_class = BBK
    # for (m,Δt) in enumerate(t_range)
    Δt=1e-3
    for int_class in [BAOAB,OBABO,BBK,GJF,VEC] # 
        params=TrajsParams(n_steps = 1e5, n_trajs = 50, n_save_iters = 1)
        println(String(Symbol(int_class))," ", Δt)
        integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
        trajs=run_trajectories(integrator; params = params,to_save=["x","v","v_mid"], init_conds_args=init_conds)

        var_meanvel=zeros(length(lambda_range),2)
        err_meanvel=zeros(length(lambda_range),2)
        var_stdvel=zeros(length(lambda_range),2)
        for (n,lambda) in enumerate(lambda_range)
            x = vcat([hcat(velocity_diff(trj.xt[1][1:end,1],trj.xt[2][1:end,1],lambda, Δt)...) for trj in trajs]...)
            var_meanvel[n,:] = mean(x,dims=1)
            cov_mat = cov(x)
            var_stdvel[n,1]=cov_mat[1,1]/(Δt)
            var_stdvel[n,2]=cov_mat[1,2]/(Δt)
            err_meanvel[n,:] =quantile(TDist(size(x)[1]-1), 1 - 0.05/2) * std(x,dims=1)/sqrt(size(x)[1])
        end
         plot!(lambda_range, var_meanvel[:,1], marker=(:circle,10),label="$(String(Symbol(int_class))) $Δt" )
         plot!(lambda_range, var_stdvel, marker=(:auto,10),label="$(String(Symbol(int_class))) std" )
    end
    xlabel!("λ")
    ylabel!("<vₜ aₜ>")
end
