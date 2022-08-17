

abstract type AbstractForce end

abstract type AbstractForceFromBasis <:AbstractForce  end

# On doit faire ensuite un type qui donne une force dérivant d'un potentiel
# Un type qui calcul des forces à partir d'une base de function

struct ForceFromPotential{F<:Function, G<: Function} <: AbstractForce
	gradV!::F
	V::G # On stocke aussi la valeur du potential
	ndim::Int
end

function ForceFromPotential(potential::String,ndim=1::Int)
	x₀=Vector{Float64}(undef,ndim)
	V = getfield(LangevinIntegrators, Symbol(potential)) # Eventuellement créer un sous module avec juste les potentials
	cfg = ForwardDiff.GradientConfig(V, x₀);
    return ForceFromPotential((gradV, x)-> ForwardDiff.gradient!(gradV, V, x, cfg), V, ndim)
end

#Se baser sur ApproxFun et (BasisMatrices.jl/ Dierckx.jl) et probablement faire un sous-type selon le package sous jacent
struct ForceFromBasis <: AbstractForceFromBasis #
	basis::Vector{Fun}
	ndim::Int
end

function ForceFromBasis(type::String,coeffs::Array{TF}) where{TF<:AbstractFloat}
	ndim=size(coeffs)[1]
	function_space=getfield(ApproxFun, Symbol(type))
	basis=Vector{Fun}(undef,ndim)
	for d in 1:ndim
		basis[d]=Fun(function_space,coeffs[d,:])
	end
	return ForceFromBasis(basis,ndim)
end

# struct ForceFromSplines <: AbstractForceFromBasis # Ca utilisera Dierckx.jl
# end

#La fonction qui calcule la force
function forceUpdate!(force::ForceFromPotential,f::Vector{TF}, x::Vector{TF})  where{TF<:AbstractFloat}
	force.gradV!(f,x)
end

function forceUpdate!(force::ForceFromBasis,f::Vector{TF}, x::Vector{TF})  where{TF<:AbstractFloat}
	for d in 1:force.ndim
		f[d]=force.basis[d](x[1]) # TODO, dealing with nd function
	end
end

#function forceUpdate!(force::ForceFromSplines,f::Vector{TF}, x::Vector{TF})  where{TF<:AbstractFloat}
#	force.gradV!(f,x)
#end


function forceUpdate(force::FP, x::Vector{TF})   where{FP<:AbstractForce, TF<:AbstractFloat}
	f=similar(x)
	forceUpdate!(force,f,x)
    return f
end
