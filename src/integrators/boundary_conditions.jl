
abstract type AbstractSpace end # A voir si on peut pas combiner avec les BC
#On aurait 2 functions apply_bc pour rester dans le manifold et
# apply_metric pour les forces de courbures

abstract type AbstractConstraints end
# Example de contraintes, quand l'état A est atteint on remet le système dans l'état B
# Différence entre BC et Contraintes, les BC s'appliquent element-wise et les contraintes globablement
# Donc en 1D pas de différence
# Les contraintes permettent de faire des BC bizarres si besoin

# On va définir l'espace sphére


#Generic function, that mean noPBC
function apply_bc(
    bc::ABC,
    x_i::TF,
) where {ABC<:AbstractBoundaryConditions,TF<:AbstractFloat}
    return xi
end


abstract type AbstractBoundaryConditions <: AbstractConstraints end

struct noBC{TF<:AbstractFloat} <: AbstractBoundaryConditions end

struct PBC{TF<:AbstractFloat} <: AbstractBoundaryConditions
    low_lim::TF
    upp_lim::TF
end

struct ReflectingBC{TF<:AbstractFloat} <: AbstractBoundaryConditions # Mais ça doit récupérer la vitesse
    low_lim::TF
    upp_lim::TF
end

#Generic function, that mean noPBC
function apply_bc(
    bc::ABC,
    x_i::TF,
) where {ABC<:AbstractBoundaryConditions,TF<:AbstractFloat}
    return xi
end

#Generic function, that mean noPBC
function apply_bc(
    bc::ABC,
    x_i::TF,
) where {ABC<:AbstractBoundaryConditions,TF<:AbstractFloat}
    return xi
end

#Il faut en fait l'appliquer element wise
function apply_bc(bc::PBC, x_i::TF) where {TF<:AbstractFloat}
    return bc.low_lim + mod(x_i - bc.low_lim, bc.upp_lim - bc.low_lim)
end
