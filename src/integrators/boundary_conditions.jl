
abstract type AbstractSpace end # A voir si on peut pas combiner avec les BC
#On aurait 2 functions apply_bc pour rester dans le manifold et
# apply_metric pour les forces de courbures

abstract type AbstractConstraints end
# Example de contraintes, quand l'état A est atteint on remet le système dans l'état B
# Différence entre BC et Contraintes, les BC s'appliquent element-wise et les contraintes globablement
# Donc en 1D pas de différence
# Les contraintes permettent de faire des BC bizarres si besoin

# On va définir l'espace sphére
function apply_space!(space::Nothing,x::Array{TF}) where {TF<:AbstractFloat}

end

function apply_space!(space::Nothing,x::Array{TF},v::Array{TF}) where {TF<:AbstractFloat}

end

"""
When we want to appy boundary conditions, define a SeparateSpace, that hold one BC or constraints per dimension.
That can hold a constraints on some dimension if needed
"""
struct SeparateSpace <: AbstractSpace
    bcs::Array{AbstractConstraints}
    dim::Int
    SeparateSpace(bcs) = new(bcs,length(bcs))
end


function apply_space!(space::SeparateSpace,x::Array{TF}) where {TF<:AbstractFloat}
    for n in 1:space.dim
        x[n], _  = apply_bc(space.bcs[n],x[n],0.0)
    end
end

function apply_space!(space::SeparateSpace,x::Array{TF},v::Array{TF}) where {TF<:AbstractFloat}
    for n in 1:space.dim
        x[n], v[n] = apply_bc(space.bcs[n],x[n],v[n])
    end
end

abstract type AbstractBoundaryConditions <: AbstractConstraints end

struct noBC <: AbstractBoundaryConditions end

function apply_bc(bc::noBC, x_i::TF, v_i::TF) where {ABC<:AbstractBoundaryConditions,TF<:AbstractFloat}
    return x_i, v_i
end

struct PBC{TF<:AbstractFloat} <: AbstractBoundaryConditions
    low_lim::TF
    upp_lim::TF
end

function apply_bc(bc::PBC, x_i::TF, v_i::TF) where {TF<:AbstractFloat}
    return bc.low_lim + mod(x_i - bc.low_lim, bc.upp_lim - bc.low_lim), v_i
end

struct ReflectingBC{TF<:AbstractFloat} <: AbstractBoundaryConditions # Mais ça doit récupérer la vitesse
    low_lim::TF
    upp_lim::TF
end

function apply_bc(bc::ReflectingBC, x_i::TF, v_i::TF) where {TF<:AbstractFloat}
    if x_i >= bc.upp_lim
        v_i=-v_i
        x_i = bc.upp_lim
    elseif x_i <= bc.low_lim
        v_i=-v_i
        x_i = bc.low_lim
    end
    return x_i,v_i
end

# For more general space, define an apply_space more general
