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
    #The various et of observer that can be defined
    obs_conf=[]
    if haskey(conf,"dump")
        dump_dict=retrieve(conf, "dump")
        dump_dict["name"]="dump"
        push!(obs_conf,dump_dict)
    end
    if haskey(conf,"logging")
        log_dict=retrieve(conf, "logging")
        log_dict["name"]="log"
        push!(obs_conf,log_dict)
    end
    if haskey(conf,"observer")
        push!(obs_conf,retrieve(conf, "observer"))
    end
    params = LangevinParams(sampling_conf,obs_conf)

    # The information about the initial conditions
    init_conds_args=Dict()
    if haskey(conf,"init_position")
        init_conds_args["position"]=retrieve(conf, "init_position")
    end
    if haskey(conf,"init_velocity")
        init_conds_args["velocity"]=retrieve(conf, "init_velocity")
    end
    if haskey(conf,"init_hidden")
        init_conds_args["hidden"]=retrieve(conf, "init_hidden")
    end
    if haskey(conf,"init_kernel")
        init_conds_args["kernel"]=retrieve(conf, "init_kernel")
    end

    #Ensuite on crée la structure qui tient la force,

    #elle qui tient l'intégrateur en fonction des paramètres
	#get(collection, key, default) pour avoir une valeur par défaut si la clé n'est pas présente
    if sampling_conf["integration"] in ["EM","euler"]
        integrator=EM(force::FP, 1.0/sampling_conf["temperature"], sampling_conf["dt"])
	end

    return integrator, params, init_conds_args
end


"""
Function to read NPZ config from GLE_analysisEM
TODO: Sortir un npz aussi de VolterraBasis et ensuite faire un read_npz_kernel
puis une function read_npz qui selectionne hidden ou gle selon le nom des key dans le npz
"""

function read_hidden_npz(file::String; integrator_type="EM"; kwargs...)
    vars = npzread("data.npz")

    #TODO some sanity checks
    dt=get(kwargs,:dt, vars["dt"]) # To allow changing the time step
    if integrator_type in ["EM","euler"] && vars["dim_h"] > 0
        integrator = EM_Hidden(force,vars["A"],vars["C"],dt,vars["dim_x"])
    elseif integrator_type == "aboba" && vars["dim_h"] > 0
        integrator = ABOBA_Hidden(force,vars["A"],vars["C"], dt,vars["dim_x"])
    end
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


struct LangevinParams where {AO <: AbstractObserver} # La dedans on stocke les trucs initialisé
    n_steps::Int
	n_trajs::Int
    #Par défaut le tableaux suivants sont vide
    observers::Array{AO}
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

function LangevinParams(sampling_dict)
    return LangevinParams(sampling_dict["n_steps"],sampling_dict["n_trajs"], [])
end

function LangevinParams(sampling_dict,obs_dict)
    obs_list=initialize_observers(obs_dict)
    return LangevinParams(sampling_dict["n_steps"],sampling_dict["n_trajs"], obs_list)
end

function addObserver(params::LangevinParams;kwargs...)
    push!(params.constraints,new_obs)
end
