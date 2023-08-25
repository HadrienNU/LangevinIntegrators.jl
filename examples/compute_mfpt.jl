#=
run_kernel.jl

Run a simple trajectory for a kernel
=#

using LangevinIntegrators
using DelimitedFiles

force = ForceFromPotential("DoubleWell",potential_kwargs=Dict(:a=>1.0))
addFix!(force,Stopper([1.0]))
params = TrajsParams(
    n_steps = 10^5,
    n_trajs = 2000,
)
Δt = 2e-3
gamma =1.0
β = 1.0
res = zeros(14,4)
res[1:end,1] = [-4.0,-3.0,-2.0,-1.5,-1.0,-0.95,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,0.95]
for (n,qval) in enumerate(res[1:end,1])
    println("------------")
    println("q ",qval)
    init_conds = Dict(
        "position" => Dict("type" => "cste", "val" => qval),
        # "velocity" => Dict("type" => "cste", "val" => sqrt(1 / β)),
        "velocity" => Dict("type" => "gaussian", "std" => sqrt(1 / β)),
    )

    integrator = VECSp(force, β,x-> gamma+5*x[1]^2, 1.0, Δt, 1)
    fpt,reached =run_fpt(integrator; params = params, init_conds_args = init_conds)

    mean_rate = sum(reached)/sum(fpt)
    # Error bound are derived from exponential distribution maximum likelihood
    println("Rate estimate: ",mean_rate," +/- ",mean_rate*exp(-1.96/sqrt(sum(reached)))," ", mean_rate*exp(1.96/sqrt(sum(reached)))) # approximate 95% CI for rate
    mfpt= 1.0/mean_rate
    println("MFPT estimate: ",mfpt," +/- ", mfpt*exp(-1.96/sqrt(sum(reached)))," ",mfpt*exp(1.96/sqrt(sum(reached)))) # approximate error bound for MFPT
    res[n,2]  = mfpt
    res[n,3]  = mfpt*exp(-1.96/sqrt(sum(reached)))
    res[n,4]  = mfpt*exp(1.96/sqrt(sum(reached)))
end
writedlm("mfpt_avg_v.dat",res)
