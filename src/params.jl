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

    return params, integrator, init_conds_args
end

function select_integrator(name::String)
    return integrator
end

"""
Function to initialize the init_cond
"""
function initialize_initcond(integrator;kwargs...)
    # En vrai ça se contente de savoir si on doit générer, 1 2 ou 3 init_cond et ça appelle get_init_conditions qui les crée
    # Ca permet de définir des valeurs par défauts si rien n'est donné
    # ça retourne un vecteur de init_cond
    intcond_pos = get_init_conditions(get(args,:position,Dict("type"=>"Cste"))
    if integrator <:OverdampedIntegrator
        # if integrator <: HiddenOverdampedIntegrator
        #     initcond_hidden = get_init_conditions(args["hidden"])
        #     return [intcond_pos,initcond_hidden]
        # elseif integrator <: KernelOverdampedIntegrator
        #     initcond_mem = get_init_conditions(args["memory"])
        #     return [intcond_pos,initcond_mem]
        # end
        return [intcond_pos]
    else
        initcond_velocity = get_init_conditions(get(args,:velocity,Dict("type"=>"Gaussian","std"=>1.0)) # à remplacer la la maxelliene
        if integrator <: HiddenIntegrator
            initcond_hidden = get_init_conditions(args["hidden"])
            return [intcond_pos,initcond_velocity,initcond_hidden]
        elseif integrator <: KernelIntegrator
            initcond_mem = get_init_conditions(args["memory"])
            return [intcond_pos,initcond_velocity,initcond_mem]
        end
        if integrator <: InertialIntegrator
            return [intcond_pos,initcond_velocity]
        end
    end
    # raise error
end

"""
Function to initialize the observers
"""
function initialize_observers(args,integrator)
    # Ca prend une liste de dict et ça génère une observable par element de la liste
    return []
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



struct LangevinParams # La dedans on stocke les trucs initialisé
    n_steps::Int
	n_trajs::Int
    observers::Array{AbstractObserver}
end

"""
    LangevinParams(;n_steps = 10^4, n_save_iters=1)
Set options for samplers.
### Fields
* n_steps       - Set the number of iterations of the sampler
* n_save_iters  - Set the frequency at which iterations are saved.  If
                  n_save_iters=1, every iteration is saved.  If n_save_iters=n_steps,
                  only the final iteration is saved.
"""
function LangevinParams(; n_steps = 10^4, n_trajs=1)
    return LangevinParams(n_steps,n_trajs, [])
end

function LangevinParams(sampling_dict,obs_dict)
    obs_list=initialize_observers(obs_dict)
    return LangevinParams(sampling_dict["n_steps"],sampling_dict["n_traj"], obs_list)
end

# function LangevinParams(sampling_dict,logging_dict,dump_dict, init_dict)
#
#     return LangevinParams(n_steps, n_save_iters, floor(Int, n_steps / n_save_iters))
# end
