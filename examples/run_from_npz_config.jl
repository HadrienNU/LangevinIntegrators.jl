#=
main.jl

Read parameters file and launch the trajectories
=#

using LangevinIntegrators

let
    npz_file= (length(ARGS) > 0 ? ARGS[1] : "coeffs.npz")
    #We should take the config file name as argument
    integrator,params,init_conf=set_hidden_from_npz(npz_file; n_steps=10^5, n_trajs=5, n_save_iters=10^4, save_filename_pattern="trajectory_*.dat")
    run_trajectories(integrator; params = params, init_conds_args=init_conf)
end
