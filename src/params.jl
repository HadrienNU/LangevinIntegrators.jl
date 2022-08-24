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
import ConfParser: retrieve


function retrieve(s::ConfParse, block::String, key::String, t::Type, default)
    if haskey(s._data[lowercase(block)],key)
        k = s._data[lowercase(block)][key]
        return parse(t, length(k) == 1 ? first(k) : k)
    else
        return default
    end
end

function check_str2(a)
    return
end


function getsubDict(s::ConfParse, block::String)
    new_dict=Dict()
    for (key,k) in s._data[lowercase(block)]
        if length(k) == 1
            new_val = tryparse(Float64, first(k)) !== nothing ?  parse(Float64,first(k)) : first(k)
        else
            new_val= k
        end
        new_dict[key] = new_val
    end
    return new_dict
end

function read_conf(file::String)
    conf = ConfParse(file)
    parse_conf!(conf)

    if !haskey(conf,"sampling")
        throw(ArgumentError("The config file should have a sampling section"))
    # elseif !haskey(conf,"physics")
    #     throw(ArgumentError("The config file should have a physics section"))
    end

    n_steps = retrieve(conf, "sampling","nsteps",Int64)
    n_trajs = retrieve(conf, "sampling","n_trajs", Int64 , 1)
    # Il faut convertir le dict qui a des array de string en Int
    #The various et of observer that can be defined
    obs_conf=[]
    if haskey(conf,"dump")
        dump_dict=getsubDict(conf, "dump")
        dump_dict["name"]="dump"
        push!(obs_conf,dump_dict)
    end
    if haskey(conf,"logging")
        log_dict=getsubDict(conf, "logging")
        log_dict["name"]="log"
        push!(obs_conf,log_dict)
    end
    if haskey(conf,"observer")
        push!(obs_conf,getsubDict(conf, "observer"))
    end
    params = LangevinParams(obs_conf;n_steps=n_steps,n_trajs=n_trajs)

    # The information about the initial conditions
    init_conds_args=Dict()
    if haskey(conf,"init_position")
        init_conds_args["position"]=getsubDict(conf, "init_position") # Es-ce qu'on ne voudrait pas convertir tous en Float
    end
    if haskey(conf,"init_velocity")
        init_conds_args["velocity"]=getsubDict(conf, "init_velocity")
    end
    if haskey(conf,"init_hidden")
        hidden_dict = getsubDict(conf, "init_hidden")
        if haskey(hidden_dict,"file") &&  split(hidden_dict["hidden"], ".")[2] =="npz" #If init_cond from .npz
            data = np.load(hidden_dict["hidden"], allow_pickle=true)
            init_conds_args["hidden"]=Dict("type"=>"Gaussian","mean"=> get(data,:µ0),"std"=>1.0)
        else
            init_conds_args["hidden"]=hidden_dict

    end
    if haskey(conf,"init_kernel")
        init_conds_args["kernel"]=getsubDict(conf, "init_kernel")
    end

    return params, init_conds_args
end

function read_integrator_conf(file::String)
    conf = ConfParse(file)
    parse_conf!(conf)

    if !haskey(conf,"sampling")
        throw(ArgumentError("The config file should have a sampling section"))
    elseif !haskey(conf,"physics")
        throw(ArgumentError("The config file should have a physics section"))
    end

    #Ensuite on crée la structure qui tient la force,
    physics_conf=getsubDict(conf, "physics")
    if haskey(physics_conf,"force") # Si il y a une clé force, elle domaine la clé potential
        #Dans ce cas, l'argument est soit un fichier, soit un tuple "name",coeffs
        if typeof(physics_conf["force"]) == String # Assume this is a file
            ext=split(physics_conf["force"], ".")[2]
            if ext =="npz"
                np = pyimport("numpy")
                data_force = np.load(file, allow_pickle=true)
                force=force_from_dict(get(data,:force),get(data,:basis)[])
            else
                throw(ArgumentError("Unable to read force file"))
            end
        else # Assume this is a type and its param
            type=lowercase(physics_conf["force"][1])
            coeffs=parse.(Float64,physics_conf["force"][2:length(physics_conf["force"])])
            if type in ["splines", "bsplines"]
                k=retrieve(conf,"physics","splines_k",Int64, 3)
                knots=parse.(Float64,physics_conf["splines_knots"])
                ForceFromSplines(k,knots,coeffs)
            else
                force=ForceFromBasis(physics_conf["force"][1],reshape(coeffs,(1,length(coeffs))))
            end
        end
    elseif haskey(physics_conf,"potential")
        if  typeof(physics_conf["potential"]) == String
            force=ForceFromPotential(physics_conf["potential"])
        else
            force=ForceFromPotential(physics_conf["potential"][1])
        end
    else
        force=nothing
    end

    #elle qui tient l'intégrateur en fonction des paramètres
	#get(collection, key, default) pour avoir une valeur par défaut si la clé n'est pas présente
    integrator_type = lowercase(retrieve(conf,"sampling","integration"))
    Δt = retrieve(conf,"sampling","dt",Float64)
    if integrator_type in ["em","euler"]
        temp= retrieve(conf,"sampling","temperature",Float64,1.0)
        integrator=EM(force, 1.0/temp, Δt)
    elseif integrator_type  in ["aboba","baoab","bbk","gjf"]
        temp= retrieve(conf,"sampling","temperature",Float64,1.0)
        γ = retrieve(conf,"physics","friction",Float64,1.0)
        mass = retrieve(conf,"physics","mass",Float64,1.0)
        integrator=getfield(LangevinIntegrators, Symbol(uppercase(integrator_type)))(force, 1.0/temp, γ, mass, Δt)
    elseif integrator_type =="verlet"
        mass = retrieve(conf,"physics","mass",Float64,1.0)
        integrator=Verlet(force, mass, Δt)
    elseif integrator_type in ["em_hidden","hidden_em","euler_hidden","hidden_euler"]
        integrator = read_integrator_hidden_npz(retrieve(conf,"physics","hidden");integrator_type="EM",force=force)
    elseif integrator_type in ["aboba_hidden","hidden_aboba"]
        integrator = read_integrator_hidden_npz(retrieve(conf,"physics","hidden");integrator_type="aboba",force=force)
	end

    return integrator
end



"""
Function to read NPZ config from GLE_analysisEM
TODO: Sortir un npz aussi de VolterraBasis et ensuite faire un read_npz_kernel
puis une function read_npz qui selectionne hidden ou gle selon le nom des key dans le npz
"""

function read_integrator_hidden_npz(file::String; integrator_type="EM", kwargs...)
    # vars = npzread(file)
    data = np.load(file, allow_pickle=true)
     # Check if this is a spline fct alors on doit passer le niveau d'après
    dim = size(get(data,:force))[1]
    dim_h = size(get(data,:A))[1]-dim
    #TODO some sanity checks, such as checking dimension
    Δt=get(kwargs,:dt, get(data,:dt)[]) # To allow changing the time step

    force_ = get(kwargs,:force,nothing)
    force = force == nothing  ? force_from_dict(get(data,:force),get(data,:basis)[]) : force
    if dim_h > 0
        if lowercase(integrator_type) in ["em","euler"]
            integrator = EM_Hidden(force,get(data,:A),get(data,:C),Δt,dim)
        elseif lowercase(integrator_type) == "aboba"
            integrator = ABOBA_Hidden(force,get(data,:A),get(data,:C), Δt,dim)
        end
    else
        integrator = ABOBA(force,get(data,:C)[1,1],get(data,:A)[1,1],1.0, Δt)
    end
    return integrator
end


function force_from_dict(coeffs,args)
    force=ForceFromPotential("Flat")
    basis_type=get(args,:basis_type)
    if basis_type == "bsplines"
        force= ForceFromScipySplines(get(args,:bsplines)[1][3],get(args,:bsplines)[1][1],coeffs)
    elseif basis_type == "linear"
        coeffs_taylor=zeros(Float64,(1,2))
        coeffs_taylor[1,:]=[get(args,:mean)[1],coeffs[1,1]]
        force=ForceFromBasis("Taylor",coeffs_taylor)
    # elseif basis_type == "polynomial"
    #
    # elseif basis_type == "bins"
    #
    elseif basis_type == "free_energy_kde"  || basis_type == "free_energy" || basis_type == "free_energy_histogram"
        force= ForceFromScipySplines(get(args,:fe_spline)[3],get(args,:fe_spline)[1],-1*coeffs[1,1].*get(args,:fe_spline)[2];der=1)
    else
        throw(ArgumentError("Unsupported basis type"))
    end
    return force
end

"""
Function to initialize the init_cond
Note si il n'y as pas ce qu'il faut ça va échouer silenciement, il faut mettre un verbose pour montrer ce qui est utilisé
"""
function initialize_initcond(integrator ,args)
    # En vrai ça se contente de savoir si on doit générer, 1 2 ou 3 init_cond et ça appelle get_init_conditions qui les crée
    # Ca permet de définir des valeurs par défauts si rien n'est donné
    # ça retourne un vecteur de init_cond
    intcond_pos = get_init_conditions(get(args,"position",Dict("type"=>"Cste")))
    if integrator isa OverdampedIntegrator
        # if integrator isa HiddenOverdampedIntegrator
        #     initcond_hidden = get_init_conditions(args["hidden"])
        #     return [intcond_pos,initcond_hidden]
        # elseif integrator isa KernelOverdampedIntegrator
        #     initcond_mem = get_init_conditions(args["memory"])
        #     return [intcond_pos,initcond_mem]
        # end
        return [intcond_pos]
    else
        initcond_velocity = get_init_conditions(get(args,"velocity",Dict("type"=>"Gaussian","std"=>1.0))) # à remplacer la la maxelliene
        if integrator isa HiddenIntegrator
            initcond_hidden = get_init_conditions(get(args,"hidden",Dict("type"=>"Gaussian","std"=>1.0)))

            return [intcond_pos,initcond_velocity,initcond_hidden]
        elseif integrator isa KernelIntegrator
            initcond_mem = get_init_conditions(get(args,"memory",Dict("type"=>"Cste")))
            return [intcond_pos,initcond_velocity,initcond_mem]
        end
        if integrator isa InertialIntegrator
            return [intcond_pos,initcond_velocity]
        end
    end
    # raise error
end



struct LangevinParams{AO <: AbstractObserver} # La dedans on stocke les trucs initialisé
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
    return LangevinParams(n_steps,n_trajs, empty([],AbstractObserver))
end


function LangevinParams(obs_dict; kwargs...)
    obs_list=initialize_observers(obs_dict)
    return LangevinParams(get(kwargs,:n_steps,10^4), get(kwargs,:n_trajs,1), obs_list)
end

function addObserver(params::LangevinParams;kwargs...)
    push!(params.constraints,new_obs)
end
