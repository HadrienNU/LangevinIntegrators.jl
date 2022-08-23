#=

Dans ce fichier on définit un ensemble de fix, de contraintes et de boundary conditions (on pourra toujours bouger les contraintes si on en a trop)

Fix: (on pourrait les définir via plumed mais comme ils sont simples ça permet d'éviter d'avoir à gérer plumed en plus)
	-UWall -> on mets des mur exponentielles pour compenser tout trop sur le fit de la force
	-LWall
	-quadratic -> on rajoute une force quadratique

	-> Le fix plumed est défini dans le module Plumed

=#


abstract type AbstractFix end

struct Wall{TF<:AbstractFloat} <:AbstractFix
	exponent::Int64 # Il faut vérifier que exponent est >1
	at::Array{TF}
	strenght::Array{TF}
	dir:: TF
end


function apply_fix!(fix::Wall,x::Array{TF},f::Array{TF}) where {TF<:AbstractFloat}
	uscaled=x.- fix.at # Dire que si c'est négatif ou positif selon dir, on met à zero
	uscaled[dir*uscaled.>0.0].=0.0 #selon la direction on met à zero les composantes positives ou négatives
	f.+= -fix.exponent*fix.strenght.*uscaled.^(fix.exponent-1)
	return 0
end


function init_fix(fix::AbstractFix; kwargs...) # Différent du constructeur au sesn que c'est appelé au début de chaque traj

end

function apply_fix!(fix::AbstractFix,x::Array{TF},f::Array{TF}) where {TF<:AbstractFloat}
	return 0
end


function close_fix(fix::AbstractFix)

end
