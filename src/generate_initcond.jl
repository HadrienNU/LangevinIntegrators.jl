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

abstract type AbstractInitCond end

abstract type AbstractRandomInitCond <: AbstractInitCond  end

struct Constant_InitCond{TF<:AbstractFloat} <: AbstractInitCond
    val::Array{TF}
    dim::Int64
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

struct Gaussian_InitCond{TF<:AbstractFloat} <: AbstractRandomInitCond
    mean::Array{TF}
    std::TF
    dim::Int64
end

# struct PMF_Init <: AbstractRandomInitCond where{TF<:AbstractFloat}
#     mean::Array{TF}
#     std::TF
#     dim::int64
# end

#Et les constructeurs de ces machins vont prendre l'intégrateur et les params en arguments, du moins un sous-ensemble des params comme ça on passe params["init_position"], params["init_velocity"]
# En vrai on va juste avoir un grand constructeur qui va s'appellait init_conditions qui retournera selon un if la bonne chose

function get_init_conditions(args::Dict)
    type=get(args,:type,"Constant")
    dim=get(args,:dim,1)
    if type in ["cste","Constante"]
        return Constant_InitCond(get(args,:val,zeros(dim)),dim)
    else
        println("Unknwon initializer")
    end
end

function generate_initcond(init_cond::Constant_InitCond;  kwargs...)
    return init_cond.val
end

function generate_initcond(init_cond::Array_InitCond; id=1, kwargs...)
    return init_cond.val[1+(id-1)%val.size[1],:]
end

function generate_initcond(init_cond::Uniform_InitCond; kwargs...)
    return low.+(high.-low).*rand(init_cond.dim)
end


function generate_initcond(init_cond::Gaussian_InitCond; kwargs...)
    return init_cond.mean .+ init_cond.std.*randn(init_cond.dim)
end
