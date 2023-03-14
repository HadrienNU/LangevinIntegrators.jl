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
function scalar_product(xt, vt, lambda,Δt, γ)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    ut0 = (xt[2:end-1]-xt[1:end-2])/Δt
    ut = lambda*(xt[3:end]-xt[2:end-1])/Δt+ (1.0 -lambda)*(xt[2:end-1]-xt[1:end-2])/Δt
    # acc = (ut[2:end]-ut[1:end-1])/ Δt
    # return ut.*(acc.+γ*ut)
    # return ut[1:end-1].*(((ut[2:end]-ut[1:end-1])/ Δt)+γ*ut[1:end-1])
    # return (((ut[2:end]-ut[1:end-1])/ Δt)+γ*ut[1:end-1]).^2
    # return (((ut[2:end]-ut[1:end-1])/ Δt)+γ*ut0[1:end-1]).*ut0[1:end-1]
    return (vt[3:end-1]-ut[2:end]).*(vt[2:end-2]-ut[1:end-1])
end


function likelihood_vec(xt, lambda,Δt, γ)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    ut1 = (xt[3:end]-xt[2:end-1])/Δt
    ut = lambda*(xt[3:end]-xt[2:end-1])/Δt+ (1.0 -lambda)*(xt[2:end-1]-xt[1:end-2])/Δt
    left = ut-ut1
    return left.*(2*ut1/Δt-acc)+γ*left.*left
    # return ut.*ut
end


function likelihood(xt, lambda,Δt, γ)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    ut = lambda*(xt[3:end]-xt[2:end-1])/Δt+ (1.0 -lambda)*(xt[2:end-1]-xt[1:end-2])/Δt
    return ((ut[2:end]+(γ-1)*ut[1:end-1]).^2)
    # return ut.*ut
end

function scalar_product_vel(xt, vt, Δt, γ)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    ut0 = (xt[2:end-1]-xt[1:end-2])/Δt
    return vt[2:end-1].*acc +γ*vt[2:end-1].*ut0
    # return vt[2:end-1].*vt[2:end-1]
end

let
    force=ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    t_range = LinRange(5e-4,5e-2,5)
    lambda_range = LinRange(-0.5,1.0,7)
    init_conds=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"gaussian", "std"=> sqrt(1/β)))
    inspectdr()

    every= 1

    plot()
    # int_class = BBK
    # for (m,Δt) in enumerate(t_range)
    Δt=1e-3
    for int_class in [BAOAB,OBABO,BBK,GJF,VEC]
        params=TrajsParams(n_steps = 1e5, n_trajs = 150, n_save_iters = 1)
        println(String(Symbol(int_class))," ", Δt)
        integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
        trajs=run_trajectories(integrator; params = params,to_save=["x","v","v_mid"], init_conds_args=init_conds)

        c2 = (1.0-integrator.c₂)/Δt

        # x = vcat([scalar_product_vel(trj.xt[1][:,1],trj.xt[2][:,1], Δt, c2) for trj in trajs]...)
        # val_velocity = mean(x)
        # err_vel = quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        # x = vcat([scalar_product_vel(trj.xt[1][:,1],trj.xt[3][:,1], Δt, c2) for trj in trajs]...)
        # val_velocity_mid = mean(x)
        # err_vel_mid = quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        # plot!(lambda_range, val_velocity*ones(length(lambda_range)), yerr =err_vel , label="Integrator velocity $(String(Symbol(int_class))) $Δt")
        # plot!(lambda_range, val_velocity_mid*ones(length(lambda_range)), yerr =err_vel_mid , label="Integrator half step velocity $(String(Symbol(int_class))) $Δt")

        var_scalar_product=zeros(size(lambda_range))
        err_scalar_product=zeros(size(lambda_range))
        for (n,lambda) in enumerate(lambda_range)
            x = vcat([scalar_product(trj.xt[1][1:every:end,1],trj.xt[2][1:every:end,1],lambda, every*Δt, c2) for trj in trajs]...)
            var_scalar_product[n] = mean(x)
            err_scalar_product[n] =quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        end
         plot!(lambda_range, var_scalar_product, yerr = err_scalar_product,marker=(:circle,1),label="$(String(Symbol(int_class))) $Δt" )
    end
    xlabel!("λ")
    ylabel!("<vₜ aₜ>")
end
