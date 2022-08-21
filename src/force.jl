

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

#Se baser sur ApproxFun et (BasisMatrices.jl/ BSplineKit.jl) et probablement faire un sous-type selon le package sous jacent
struct ForceFromBasis <: AbstractForceFromBasis #
	basis::Vector{Fun}
	ndim::Int
end

function ForceFromBasis(type::String,coeffs::Array{TF}) where{TF<:AbstractFloat}
	if ndims(coeffs) >=2
		ndim=size(coeffs)[1]
	else
		ndim=1
	end
	function_space=getfield(ApproxFun, Symbol(type))
	basis=Vector{Fun}(undef,ndim)
	for d in 1:ndim
		basis[d]=Fun(function_space,coeffs[d,:])
	end
	return ForceFromBasis(basis,ndim)
end

# The struct is writed to have the same interface than Force from basis to have only one forceUpdate
struct ForceFromSplines <: AbstractForceFromBasis # use BSplineKit
	basis::Array#{Splines}
	ndim::Int
end

function ForceFromSplines(k::Int,knots::AbstractVector,coeffs::Array{TF}) where{TF<:AbstractFloat}
	if ndims(coeffs) >=2
		ndim=size(coeffs)[1]
	else
		ndim=1
	end
	#Note the order of the splines is k+1 with k the degree
	B = BSplineBasis(BSplineOrder(k+1), knots;augment=false);
	basis=Vector{Splines}(undef,ndim)
	for d in 1:ndim
		basis[d]=Spline(B, coeffs[d,:])
	end
	return ForceFromSplines(basis,ndim)
end

#La fonction qui calcule la force a actualiser pour n'avoir que forceUpdate! from basis
function forceUpdate!(force::ForceFromPotential,f::Vector{TF}, x::Vector{TF})  where{TF<:AbstractFloat}
	force.gradV!(f,x)
	f.*=-1
end

function forceUpdate!(force::FB,f::Vector{TF}, x::Vector{TF})  where{FB<:AbstractForceFromBasis,TF<:AbstractFloat}
	for d in 1:force.ndim
		f[d]=force.basis[d](x[1]) # TODO, dealing with nd function
	end
end

function forceUpdate(force::FP, x::Vector{TF})   where{FP<:AbstractForce, TF<:AbstractFloat}
	f=similar(x)
	forceUpdate!(force,f,x)
    return f
end

# A function to plot value of the force for comparison with python code
function force_eval(force::FP, x::Vector{TF}) where{FP<:AbstractForce,TF<:AbstractFloat}
	f=similar(x)
	for (n, val) in enumerate(x)
		f[n] = forceUpdate(force,[val])[1]
	end
	return f
end
