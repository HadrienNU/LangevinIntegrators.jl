#=
main.jl

Compute diffusion coefficients for various set of integrators and compare it to expected value
=#

using Plots
using LangevinIntegrators

let
    force=ForceFromPotential("Flat")
    γ = 1.0
    β = 1.0
    length_trj = 100.0
    t_range = LinRange(1e-4,1e-1,5)
    inspectdr()

    for int_class in [ABOBA,BAOAB,BBK,GJF]

        diff_coeff=zeros(size(t_range))
        for (n,Δt) in enumerate(t_range)
            params=TrajsParams(n_steps = floor(length_trj/Δt), n_trajs = 250, n_save_iters = floor(length_trj/Δt))
            println(String(Symbol(int_class))," ", Δt)

            integrator=int_class(force, β , γ, 1.0, Δt, 1)
            trajs=run_trajectories(integrator; params = params)
            for trj in trajs
                diff_coeff[n] +=(trj.xt[end])^2/(2*params.n_steps*Δt*params.n_trajs) # All trajectories start from 0
            end
        end
         plot!(t_range, diff_coeff,label=String(Symbol(int_class)))
    end
    xlabel!("Δt")
    ylabel!("<(x_t-x_0)^2>")
end
