#=

Dans ce fichier on définit un ensemble de fix, de contraintes et de boundary conditions (on pourra toujours bouger les contraintes si on en a trop)

Fix: (on pourrait les définir via plumed mais comme ils sont simples ça permet d'éviter d'avoir à gérer plumed en plus)
	-UWall -> on mets des mur exponentielles pour compenser tout trop sur le fit de la force
	-LWall
	-quadratic -> on rajoute une force quadratique

	-> Le fix plumed est défini dans le module Plumed

=#


abstract type AbstractFix end

struct Wall{TF<:AbstractFloat} <: AbstractFix
    exponent::Int64 # Il faut vérifier que exponent est >1
    at::Array{TF}
    strenght::Array{TF}
    dir::TF
    bias_energy::Float64
end

function UWall(exponent, at, strenght)
    return Wall(exponent, at, strenght, -1.0, 0.0)
end

function LWall(exponent, at, strenght)
    return Wall(exponent, at, strenght, 1.0, 0.0)
end

function apply_fix!(fix::Wall, x::Array{TF}, f::Array{TF}) where {TF<:AbstractFloat}
    uscaled = x .- fix.at # Dire que si c'est négatif ou positif selon dir, on met à zero
    uscaled[fix.dir*uscaled.>0.0] .= 0.0 #selon la direction on met à zero les composantes positives ou négatives
    f .+= -fix.exponent * fix.strenght .* uscaled .^ (fix.exponent - 1)
    return false
end

struct PolynomialForce{TF<:AbstractFloat} <: AbstractFix
    exponent::Int64 # Il faut vérifier que exponent est >1
    at::Array{TF}
    strenght::Array{TF}
    bias_energy::Float64
end


function Quadratic(at::Array{TF}, strenght::Array{TF}) where {TF<:AbstractFloat}
    return PolynomialForce(2, at::Array{TF}, strenght::Array{TF}, 0.0)
end

function apply_fix!(fix::PolynomialForce, x::Array{TF}, f::Array{TF}) where {TF<:AbstractFloat}
    f .+= -fix.exponent * fix.strenght .* (x .- fix.at) .^ (fix.exponent - 1)
    return 0
end


function init_fix(fix::AbstractFix; kwargs...) # Différent du constructeur au sesn que c'est appelé au début de chaque traj
    return fix
end

#Write a getter for the extra energy of the fixes, so when the force is asked to provide the energy it can iterate over the fixes and compute the nergy

function close_fix(fix::AbstractFix)

end
