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
function scalar_product(xt, lambda,Δt, γ)
    acc = (xt[3:end]-2*xt[2:end-1]+xt[1:end-2])/ (Δt^2)
    ut0 = (xt[2:end-1]-xt[1:end-2])/Δt
    ut = lambda*(xt[3:end]-xt[2:end-1])/Δt+ (1.0 -lambda)*(xt[2:end-1]-xt[1:end-2])/Δt
    return ut.*acc+γ*ut.*ut0
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
    cg_range = range(1,10)
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

        c2 = 0.0 # (1.0-integrator.c₂)/Δt

        val_velocity=zeros(size(cg_range))
        err_vel=zeros(size(cg_range))

        val_velocity_mid=zeros(size(cg_range))
        err_vel_mid=zeros(size(cg_range))

        var_scalar_product=zeros(size(cg_range))
        err_scalar_product=zeros(size(cg_range))
        for (n,cg_step) in enumerate(cg_range)

            x = vcat([scalar_product_vel(trj.xt[1][1:cg_step:end,1],trj.xt[2][1:cg_step:end,1], cg_step*Δt, c2) for trj in trajs]...)
            val_velocity[n] = mean(x)
            err_vel[n] = quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
            x = vcat([scalar_product_vel(trj.xt[1][1:cg_step:end,1],trj.xt[3][1:cg_step:end,1], cg_step*Δt, c2) for trj in trajs]...)
            val_velocity_mid[n] = mean(x)
            err_vel_mid[n] = quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))


            x = vcat([scalar_product(trj.xt[1][1:cg_step:end,1], 0.0, Δt*cg_step, c2) for trj in trajs]...)
            var_scalar_product[n] = mean(x)
            err_scalar_product[n] =quantile(TDist(length(x)-1), 1 - 0.05/2) * std(x)/sqrt(length(x))
        end
        plot!(cg_range, val_velocity, yerr =err_vel , label="Integrator velocity $(String(Symbol(int_class))) $Δt")
        plot!(cg_range, val_velocity_mid, yerr =err_vel_mid , label="Integrator half step velocity $(String(Symbol(int_class))) $Δt")
        plot!(cg_range, var_scalar_product, yerr = err_scalar_product,marker=(:circle,3),label="$(String(Symbol(int_class))) $Δt" )
    end
    xlabel!("n")
    ylabel!("<vₜ aₜ>")
end
