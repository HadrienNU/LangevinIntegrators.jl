#=
The struct that hold global parameters
    Le reste des parametres est passé à l'intégrateur

    Parametres:
    -Globaux
        Nombre de trajectoires
        Nombre de pas de temps
        Burning steps (mais ça va avec le logging)
        Dimension du système
        PBC (boite)
    -Logging:
        Fréquence de sortie des infos
        Si on affiche quelque chose sur l'écran

    ---> params

    -Initialisation:
        Position
        Vitesse
        Variables cachées
        Mémoire

    --> initializer

    -Intégrateur:
        Pas de temps
        Force
        Paramètres friction
        Paramètres bain

    --> integrateur

    Les fonctions, read_conf doivent donc returner un params, un intégrateur et un initializer
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

    return params, integrator
end

function select_integrator(name::String)
    return integrator
end

"""
Function to read NPZ config from GLE_analysisEM
"""

function read_npz(file::String; integrator_type="EM", dt=1.0)
    vars = npzread("data.npz")

    #TODO some sanity checks and deal with the possibility of changing the time step

    if integrator_type in ["EM","euler"] && vars["dim_h"] > 0
        integrator = EM_Hidden(force,vars["A"],vars["C"], vars["dt"],vars["dim_x"])
    elseif integrator_type == "aboba" && vars["dim_h"] > 0
        integrator = ABOBA_Hidden(force,vars["A"],vars["C"], vars["dt"],vars["dim_x"])
    end
    return integrator
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
function LangevinParams(; n_iters = 10^4,n_trajs=1, n_save_iters = 1)

    return LangevinParams(n_iters,n_trajs, n_save_iters, floor(Int, n_iters / n_save_iters))
end

# function LangevinParams(sampling_dict,logging_dict,dump_dict, init_dict)
#
#     return LangevinParams(n_iters, n_save_iters, floor(Int, n_iters / n_save_iters))
# end
