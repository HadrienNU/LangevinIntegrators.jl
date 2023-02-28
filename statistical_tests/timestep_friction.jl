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
function scalar_product(xt, lambda,Δt)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    ut = lambda*(xt[3:end]-xt[2:end-1])/Δt+ (1.0 -lambda)*(xt[2:end-1]-xt[1:end-2])/Δt
    return ut.*acc
    # return ut.*ut
end

function scalar_product_vel(xt, vt, Δt)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    return vt[2:end-1].*acc
    # return vt[2:end-1].*vt[2:end-1]
end

let
    force=ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    t_range = LinRange(5e-4,5e-2,5)
    lambda_range = LinRange(0.0,1.0,5)
    init_conds=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"gaussian", "std"=> sqrt(1/β)))
    inspectdr()



    plot()
    # int_class = BBK
    # for (m,Δt) in enumerate(t_range)
    Δt=1e-3
    for int_class in [BAOAB,OBABO,BBK,GJF,VEC]
        params=TrajsParams(n_steps = 1e5, n_trajs = 150, n_save_iters = 1)
        println(String(Symbol(int_class))," ", Δt)
        integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
        trajs=run_trajectories(integrator; params = params,to_save=["x","v","v_mid"], init_conds_args=init_conds)

        x = vcat([scalar_product_vel(trj.xt[1][:,1],trj.xt[2][:,1], Δt) for trj in trajs]...)
        val_velocity = mean(x)
        err_vel = quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        x = vcat([scalar_product_vel(trj.xt[1][:,1],trj.xt[3][:,1], Δt) for trj in trajs]...)
        val_velocity_mid = mean(x)
        err_vel_mid = quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        plot!(lambda_range, val_velocity*ones(length(lambda_range)), yerr =err_vel , label="Integrator velocity $(String(Symbol(int_class))) $Δt")
        plot!(lambda_range, val_velocity_mid*ones(length(lambda_range)), yerr =err_vel_mid , label="Integrator half step velocity $(String(Symbol(int_class))) $Δt")

        var_scalar_product=zeros(size(lambda_range))
        err_scalar_product=zeros(size(lambda_range))
        for (n,lambda) in enumerate(lambda_range)
            x = vcat([scalar_product(trj.xt[1][:,1],lambda, Δt) for trj in trajs]...)
            var_scalar_product[n] = mean(x)
            err_scalar_product[n] =quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        end
         plot!(lambda_range, var_scalar_product, yerr = err_scalar_product,marker=(:circle,1),label="$(String(Symbol(int_class))) $Δt" )
    end
    xlabel!("λ")
    ylabel!("<vₜ aₜ>")
end
