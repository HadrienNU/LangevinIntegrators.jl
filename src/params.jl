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
    if haskey(s._data[lowercase(block)], key)
        k = s._data[lowercase(block)][key]
        return parse(t, length(k) == 1 ? first(k) : k)
    else
        return default
    end
end


function getsubDict(s::ConfParse, block::String)
    new_dict = Dict()
    for (key, k) in s._data[lowercase(block)]
        if length(k) == 1
            new_val =
                tryparse(Float64, first(k)) !== nothing ? parse(Float64, first(k)) :
                first(k)
        else
            new_val = k
        end
        new_dict[key] = new_val
    end
    return new_dict
end

function transform_dict_into_kwargs(dict)
    return Dict(Symbol(k) => v for (k, v) in dict)
end

function read_conf(file::String)
    conf = ConfParse(file)
    parse_conf!(conf)

    if !haskey(conf, "sampling")
        throw(ArgumentError("The config file should have a sampling section"))
        # elseif !haskey(conf,"physics")
        #     throw(ArgumentError("The config file should have a physics section"))
    end

    n_steps = retrieve(conf, "sampling", "nsteps", Int64)
    n_trajs = retrieve(conf, "sampling", "n_trajs", Int64, 1)
    verbose = retrieve(conf, "sampling", "verbose", Int64, 0)

    # Some informations to save the trajectories
    dump_dict = haskey(conf, "dump") ? getsubDict(conf, "dump") : Dict()
    params = TrajsParams(;
        n_steps = n_steps,
        n_trajs = n_trajs,
        verbose = verbose,
        transform_dict_into_kwargs(dump_dict)...,
    )

    # The information about the initial conditions
    init_conds_args = Dict()
    if haskey(conf, "init_position")
        init_conds_args["position"] = getsubDict(conf, "init_position") # Es-ce qu'on ne voudrait pas convertir tous en Float
    end
    if haskey(conf, "init_velocity")
        init_conds_args["velocity"] = getsubDict(conf, "init_velocity")
    end
    if haskey(conf, "init_hidden")
        hidden_dict = getsubDict(conf, "init_hidden")
        if haskey(hidden_dict, "file") && split(hidden_dict["file"], ".")[2] == "npz" #If init_cond from .npz
            data = np.load(hidden_dict["file"], allow_pickle = true)
            init_conds_args["hidden"] =
                Dict("type" => "Gaussian", "mean" => get(data, "µ_0"), "std" => 1.0)
        else
            init_conds_args["hidden"] = hidden_dict
        end

    end
    if haskey(conf, "init_kernel")
        init_conds_args["kernel"] = getsubDict(conf, "init_kernel")
    end

    return params, init_conds_args
end

function read_integrator_conf(file::String)
    conf = ConfParse(file)
    parse_conf!(conf)

    if !haskey(conf, "sampling")
        throw(ArgumentError("The config file should have a sampling section"))
    elseif !haskey(conf, "physics")
        throw(ArgumentError("The config file should have a physics section"))
    elseif !haskey(conf, "space")
        throw(ArgumentError("The config file should have a space section"))
    end

    #Dimensionality information on the system
    dim = retrieve(conf, "space", "ndim", Int64)

    #Ensuite on crée la structure qui tient la force,
    physics_conf = getsubDict(conf, "physics")
    if haskey(physics_conf, "force") # Si il y a une clé force, elle domaine la clé potential
        #Dans ce cas, l'argument est soit un fichier, soit un tuple "name",coeffs
        if typeof(physics_conf["force"]) == String # Assume this is a file
            ext = split(physics_conf["force"], ".")[2]
            if ext == "npz"
                np = pyimport("numpy")
                data_force = np.load(file, allow_pickle = true)
                force = force_from_dict(get(data_force, :force), get(data_force, :basis)[])
            else
                throw(ArgumentError("Unable to read force file"))
            end
        else # Assume this is a type and its param
            type = lowercase(physics_conf["force"][1])
            coeffs = parse.(Float64, physics_conf["force"][2:length(physics_conf["force"])])
            if type in ["splines", "bsplines"]
                k = retrieve(conf, "physics", "splines_k", Int64, 3)
                knots = parse.(Float64, physics_conf["splines_knots"])
                ForceFromSplines(k, knots, coeffs)
            else
                force = ForceFromBasis(
                    physics_conf["force"][1],
                    reshape(coeffs, (1, length(coeffs))),
                )
            end
        end
    elseif haskey(physics_conf, "potential")
        if typeof(physics_conf["potential"]) == String
            force = ForceFromPotential(physics_conf["potential"], dim)
        else
            force = ForceFromPotential(physics_conf["potential"][1], dim)
        end
    else
        force = nothing
    end

    #elle qui tient l'intégrateur en fonction des paramètres
    #get(collection, key, default) pour avoir une valeur par défaut si la clé n'est pas présente
    integrator_type = lowercase(retrieve(conf, "sampling", "integration"))
    Δt = retrieve(conf, "sampling", "dt", Float64)
    if integrator_type in ["em", "euler"]
        temp = retrieve(conf, "sampling", "temperature", Float64, 1.0)
        integrator = EM(force, 1.0 / temp, Δt, dim)
    elseif integrator_type in ["aboba", "baoab", "bbk", "gjf"]
        temp = retrieve(conf, "physics", "temperature", Float64, 1.0)
        γ = retrieve(conf, "physics", "friction", Float64, 1.0)
        mass = retrieve(conf, "physics", "mass", Float64, 1.0)
        integrator = getfield(LangevinIntegrators, Symbol(uppercase(integrator_type)))(
            force,
            1.0 / temp,
            γ,
            mass,
            Δt,
            dim,
        )
    elseif integrator_type == "verlet"
        mass = retrieve(conf, "physics", "mass", Float64, 1.0)
        integrator = Verlet(force, mass, Δt, dim)
    elseif integrator_type in ["em_hidden", "hidden_em", "euler_hidden", "hidden_euler"]
        integrator = read_integrator_hidden_npz(
            retrieve(conf, "physics", "hidden");
            integrator_type = "EM",
            force = force,
        )
    elseif integrator_type in ["aboba_hidden", "hidden_aboba"]
        integrator = read_integrator_hidden_npz(
            retrieve(conf, "physics", "hidden");
            integrator_type = "aboba",
            force = force,
        )
    end

    return integrator
end


function set_from_conf(file::String)
    params, init_conf = read_conf(file)
    integrator = read_integrator_conf(file)
    return integrator, params, init_conf
end


"""
Function to read NPZ config from GLE_analysisEM
TODO: Sortir un npz aussi de VolterraBasis et ensuite faire un read_npz_kernel
puis une function read_npz qui selectionne hidden ou gle selon le nom des key dans le npz
"""

function set_hidden_from_npz(file::String; kwargs...)
    integrator = read_integrator_hidden_npz(file; kwargs...)

    params = TrajsParams(; kwargs...)

    # The information about the initial conditions
    init_conds_args = Dict()
    init_conds_args["position"] = get(kwargs, :init_position, Dict("type" => "Cste"))
    init_conds_args["velocity"] =
        get(kwargs, :init_velocity, Dict("type" => "Gaussian", "mean" => 0.0, "std" => 1.0))
    data = np.load(file, allow_pickle = true)
    init_conds_args["hidden"] =
        Dict("type" => "Gaussian", "mean" => get(data, "µ_0"), "std" => 1.0)

    return integrator, params, init_conds_args
end


function read_integrator_hidden_npz(file::String; integrator_type = "EM", kwargs...)
    # vars = npzread(file)
    data = np.load(file, allow_pickle = true)
    # Check if this is a spline fct alors on doit passer le niveau d'après
    dim = size(get(data, :force))[1]
    dim_h = size(get(data, :A))[1] - dim
    #TODO some sanity checks, such as checking dimension
    Δt = get(kwargs, :dt, get(data, :dt)[]) # To allow changing the time step

    force = get(kwargs, :force, nothing)
    force =
        isnothing(force) ? force_from_dict(get(data, :force), get(data, :basis)[]) : force
    if dim_h > 0
        if lowercase(integrator_type) in ["em", "euler"]
            integrator = EM_Hidden(force, get(data, :A), get(data, :C), Δt, dim)
        elseif lowercase(integrator_type) == "aboba"
            integrator = ABOBA_Hidden(force, get(data, :A), get(data, :C), Δt, dim)
        end
    else
        integrator = ABOBA(force, get(data, :C)[1, 1], get(data, :A)[1, 1], 1.0, Δt)
    end
    return integrator
end


function force_from_dict(coeffs, args)
    force = ForceFromPotential("Flat")
    basis_type = get(args, :basis_type)
    if basis_type == "bsplines"
        force =
            ForceFromSplines(get(args, :bsplines)[1][3], get(args, :bsplines)[1][1], coeffs)
    elseif basis_type == "linear"
        coeffs_taylor = zeros(Float64, (1, 2))
        coeffs_taylor[1, :] = [get(args, :mean)[1], coeffs[1, 1]]
        force = ForceFromBasis("Taylor", coeffs_taylor)
        # elseif basis_type == "polynomial"
        #
        # elseif basis_type == "bins"
        #
    elseif basis_type == "free_energy_kde" ||
           basis_type == "free_energy" ||
           basis_type == "free_energy_histogram"
        force = ForceFromSplines(
            get(args, :fe_spline)[3],
            get(args, :fe_spline)[1],
            -1 * coeffs[1, 1] .* get(args, :fe_spline)[2];
            der = 1,
        )
    else
        throw(ArgumentError("Unsupported basis type"))
    end
    return force
end


struct TrajsParams # La dedans on stocke les trucs initialisé
    n_steps::Int
    n_trajs::Int
    # Some infos to save the trajectories
    n_save_iters::Int
    save_filename_pattern::Union{String,Nothing}

    verbose::Int
end

"""
    TrajsParams(;n_steps = 10^4, n_save_iters=1)
Set options for samplers.
### Fields
* n_steps       - Set the number of iterations of the integrator
* n_trajs       - Set the number of trajectories to run
* n_save_iters  - Set the frequency at which iterations are saved.  If
                  n_save_iters=1, every iteration is saved.  If n_save_iters=n_steps,
                  only the final iteration is saved.
* save_filename_pattern - The pattern for saving the trajectories. If this include a *, that will be replaced by the id of the traj.
* verbose       - Verbose argument
"""
function TrajsParams(;
    n_steps = 10^4,
    n_trajs = 1,
    n_save_iters = nothing,
    save_filename_pattern = nothing,
    verbose = 0,
    kwargs...,
)
    return TrajsParams(
        n_steps,
        n_trajs,
        isnothing(n_save_iters) ? n_steps : n_save_iters,
        save_filename_pattern,
        verbose,
    )
end
