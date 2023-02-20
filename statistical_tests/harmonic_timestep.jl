#=
main.jl

Compute diffusion coefficients for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators

let
    force=ForceFromPotential("Harmonic")
    γ = 1.0
    β = 1.0
    length_trj = 100.0
    t_range = LinRange(1e-3,1e-1,5)
    # inspectdr()
    plot(t_range, ones(length(t_range)),xaxis=:log, label="Theory")
    for int_class in [ABOBA,BAOAB,BBK,GJF]

        var_position=zeros(size(t_range))
        err_var=zeros(size(t_range))
        for (n,Δt) in enumerate(t_range)
            params=TrajsParams(n_steps = floor(length_trj/Δt), n_trajs = 250, n_save_iters = floor(0.001*length_trj/Δt))
            println(String(Symbol(int_class))," ", Δt)

            integrator=int_class(force, β , γ, 1.0, Δt, 1) # Change also initial conidition
            trajs=run_trajectories(integrator; params = params)
            for trj in trajs
                var_position[n] +=sum(trj.xt[1][:,1].^2)/(length(trj.xt[1][:,1])*params.n_trajs-1)
            end
        end
         plot!(t_range, var_position,xaxis=:log,marker=(:circle,1),label=String(Symbol(int_class)))
    end
    xlabel!("Δt")
    ylabel!("<x_t^2>")
end
