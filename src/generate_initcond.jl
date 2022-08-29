#=
initialization.jl

To get initial state for the trajectories

Il faut créer une structure qui contient les paramètres de l'initialisation

Type d'initialisation:
    -Constant
    - Array (on fournit un tableau et ça génère successivement les CI)
    - Random:
        Uniform
        Gaussien
        Selon une pmf (si on a un potentiel et des bornes, dans ce cas ça dérive de Uniform)

En fait on a 2 types d'initialisation en série et en parallèle.
Ca veut dire qu'il y a un état interne qui détermine l'initialisation d'après
Ca sera une seed pour le cas uniforme et un index pour le tableau
Pour se simplifier la vie, on peut supposer qu'il n'y a que des initializer serial et ce sont des états qu'on passe au cas distribués.
Ca permet un meileur controle de ce qu'on passe, on passe un id qui permet de gérer les 2 cas


=#

"""
Function to initialize the init_cond
Note si il n'y as pas ce qu'il faut ça va échouer silenciement, il faut mettre un verbose pour montrer ce qui est utilisé
"""
function initialize_initcond(integrator, args; verbose = 0)
    # En vrai ça se contente de savoir si on doit générer, 1 2 ou 3 init_cond et ça appelle get_init_conditions qui les crée
    # Ca permet de définir des valeurs par défauts si rien n'est donné
    intcond_pos = get_init_conditions(get(args, "position", Dict("type" => "Cste")))
    if integrator isa OverdampedIntegrator
        if integrator isa HiddenOverdampedIntegrator
            initcond_hidden = get_init_conditions(args["hidden"])
            return [intcond_pos, initcond_hidden]
        elseif integrator isa KernelOverdampedIntegrator
            initcond_mem = get_init_conditions(args["memory"])
            return [intcond_pos, initcond_mem]
        end
        return [intcond_pos]
    else
        initcond_velocity = get_init_conditions(get(args, "velocity", Dict("type" => "Gaussian", "std" => 1.0))) # à remplacer la la maxelliene
        if integrator isa HiddenIntegrator
            initcond_hidden = get_init_conditions(get(args, "hidden", Dict("type" => "Gaussian", "std" => 1.0)))

            return [intcond_pos, initcond_velocity, initcond_hidden]
        elseif integrator isa KernelIntegrator
            initcond_mem = get_init_conditions(get(args, "memory", Dict("type" => "Cste")))
            return [intcond_pos, initcond_velocity, initcond_mem]
        end
        if integrator isa InertialIntegrator
            return [intcond_pos, initcond_velocity]
        end
    end
    # raise error
    throw(ArgumentError("The integrator does not allow for initial conditons"))
end


abstract type AbstractInitCond end

abstract type AbstractRandomInitCond <: AbstractInitCond end

struct Constant_InitCond{TF<:AbstractFloat} <: AbstractInitCond
    val::Array{TF}
    dim::Int64
end

function Constant_InitCond(val::TF) where {TF<:AbstractFloat}
    return Constant_InitCond([val], 1)
end

function Constant_InitCond(val::TF, dim) where {TF<:AbstractFloat}
    return Constant_InitCond([val], 1)
end

struct Array_InitCond{TF<:AbstractFloat} <: AbstractInitCond
    val::Matrix{TF}
    dim::Int64
end

struct Uniform_InitCond{TF<:AbstractFloat} <: AbstractRandomInitCond
    low::Array{TF}
    high::Array{TF}
    dim::Int64
end

function Uniform_InitCond(low::TF, high::TF) where {TF<:AbstractFloat}
    return Uniform_InitCond([low], [high], 1)
end

function Uniform_InitCond(low::TF, high::TF, dim) where {TF<:AbstractFloat}
    return Uniform_InitCond([low], [high], 1)
end

struct Gaussian_InitCond{TF<:AbstractFloat} <: AbstractRandomInitCond
    mean::Array{TF}
    std::TF
    dim::Int64
end

function Uniform_InitCond(mean::TF, std::TF) where {TF<:AbstractFloat}
    return Uniform_InitCond([mean], [high], 1)
end

function Gaussian_InitCond(mean::TF, std::TF, dim) where {TF<:AbstractFloat}
    return Gaussian_InitCond([mean], [std], 1)
end

# struct PMF_Init <: AbstractRandomInitCond where{TF<:AbstractFloat}
#     mean::Array{TF}
#     std::TF
#     dim::int64
# end

#Et les constructeurs de ces machins vont prendre l'intégrateur et les params en arguments, du moins un sous-ensemble des params comme ça on passe params["init_position"], params["init_velocity"]
# En vrai on va juste avoir un grand constructeur qui va s'appellait init_conditions qui retournera selon un if la bonne chose

function get_init_conditions(args::Dict)
    type = get(args, "type", "Constant")
    dim = get(args, "dim", 1)
    if lowercase(type) in ["cste", "constant"]
        return Constant_InitCond(get(args, "val", zeros(dim)), dim)
    elseif lowercase(type) in ["uniform", "random"]
        low = get(args, "low", zeros(dim))
        if (low isa Number)
            low = low * ones(dim)
        else
            low = parse.(Float64, low)
        end
        high = get(args, "high", ones(dim))
        if (high isa Number)
            high = high * ones(dim)
        else
            high = parse.(Float64, high)
        end
        return Uniform_InitCond(low, high, dim)
    elseif lowercase(type) in ["gaussian", "normal"]
        mean = get(args, "mean", zeros(dim))
        if (mean isa Number)
            mean = mean * ones(dim)
        elseif mean[1] isa String
            mean = parse.(Float64, mean)
        end
        std = get(args, "std", 1.0)
        return Gaussian_InitCond(mean, std, dim)
    else
        @warn "Unknwon initializer"
    end
end

function generate_initcond(init_cond::Constant_InitCond; kwargs...)
    return init_cond.val
end

function generate_initcond(init_cond::Array_InitCond; id = 1, kwargs...)
    return init_cond.val[1+(id-1)%val.size[1], :]
end

function generate_initcond(init_cond::Uniform_InitCond; kwargs...)
    return low .+ (high .- low) .* rand(init_cond.dim)
end


function generate_initcond(init_cond::Gaussian_InitCond; kwargs...)
    return init_cond.mean .+ init_cond.std .* randn(init_cond.dim)
end
