#=
The struct that hold global parameters
    Le reste des parametres est passé à l'intégrateur
=#

function read_conf(file::String)

    conf = ConfParse("confs/config.ini")
    parse_conf!(conf)

    sampling_conf  = retrieve(conf, "sampling")

    params = LangevinParams(sampling_conf, retrieve(conf, "logging"), retrieve(conf, "dump"),retrieve(conf, "init"))

    #Ensuite on crée la structure qui tient la force,

    #elle qui tient l'intégrateur en fonction des paramètres
	#get(collection, key, default) pour avoir une valeur par défaut si la clé n'est pas présente
    if sampling_conf["integration"] in ["EM","euler"]
        integrator=EM(force::FP, 1.0/sampling_conf["temperature"], sampling_conf["dt"])
	end

    return params,integrator
end


struct LangevinParams
    n_iters::Int
	n_trajs::Int
    n_save_iters::Int
    n_save::Int
end

"""
    LangevinParams(;n_iters = 10^4, n_save_iters=1)
Set options for samplers.
### Fields
* n_iters       - Set the number of iterations of the sampler
* n_save_iters  - Set the frequency at which iterations are saved.  If
                  n_save_iters=1, every iteration is saved.  If n_save_iters=n_iters,
                  only the final iteration is saved.
"""
function LangevinParams(; n_iters = 10^4, n_save_iters = 1)

    return LangevinParams(n_iters,1, n_save_iters, floor(Int, n_iters / n_save_iters))
end

# function LangevinParams(sampling_dict,logging_dict,dump_dict, init_dict)
#
#     return LangevinParams(n_iters, n_save_iters, floor(Int, n_iters / n_save_iters))
# end
