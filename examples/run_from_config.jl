#=
main.jl

Read parameters file and launch the trajectories
=#

using LangevinIntegrators
# using ArgParse?? A bit overkilling to juste take a filename as argument

let
    #We should take the config file name as argument
    integrator,params,init_conf=read_conf("onetraj.ini")
    run_trajectories(integrator; params = params, init_conds_args=init_conf)
end
