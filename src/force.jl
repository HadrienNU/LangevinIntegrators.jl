
#S'il y a un object plumed il doit être stocké dans la force?
# En tout cas, force update doit appeler plumed_step
# La force doit stocker une liste de Fix, dont celui de plumed
# Function a créer

abstract type AbstractForce end

abstract type AbstractForceFromBasis <:AbstractForce  end

# On doit faire ensuite un type qui donne une force dérivant d'un potentiel
# Un type qui calcul des forces à partir d'une base de function

struct ForceFromPotential{F<:Function, G<: Function} <: AbstractForce
	gradV!::F
	V::G # On stocke aussi la valeur du potential
	ndim::Int
	# fixes::Array{AFix} # Uneliste de fix à appliquer
	# plumed_object  à part?
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

function ForceFromSplines(k::Int,knots::Array{TF},coeffs::Array{TF}) where{TF<:AbstractFloat}
	if ndims(coeffs) >=2
		ndim=size(coeffs)[1]
		nb_coeffs=size(coeffs)[2]
	else
		ndim=1
		nb_coeffs=length(coeffs)
		coeffs = reshape(coeffs,(1,nb_coeffs))
	end
	coeffs = coeffs[1: nb_coeffs==length(knots) ? nb_coeffs-(k+1) : nb_coeffs]
	#Note the order of the splines is k+1 with k the degree
	B = BSplineBasis(BSplineOrder(k+1), knots;augment=Val(false));
	basis=Vector{Spline}(undef,ndim)
	for d in 1:ndim
		basis[d]=Spline(B, coeffs[d,:])
	end
	return ForceFromSplines(basis,ndim)
end
#
# struct ForceFromScipySplines <: AbstractForceFromBasis  # use Scipy
# 	basis::Array#{Splines}
# 	ndim::Int
# end
#
# function ForceFromScipySplines(k::Int,knots::AbstractVector,coeffs::Array{TF}) where{TF<:AbstractFloat}
# 	if ndims(coeffs) >=2
# 		ndim=size(coeffs)[1]
# 	else
# 		ndim=1
# 	end
# 	basis=Vector{Function}(undef,ndim)
# 	for d in 1:ndim
# 		basis[d]= x-> scipy_interpolate.splev(x, (knots, coeffs[d,:], k))[]
# 	end
# 	return ForceFromScipySplines(basis,ndim)
# end

# Another version
struct ForceFromScipySplines{TF<:AbstractFloat} <: AbstractForceFromBasis  # use Scipy
	knots::Array{TF}
	coeffs_splines::Array{TF}
	k::Int
	der::Int
	ndim::Int
end
function ForceFromScipySplines(k::Int,knots::Array{TF},coeffs::Array{TF}; der=0)  where{TF<:AbstractFloat}
	if ndims(coeffs) >=2
		ndim=size(coeffs)[1]
	else
		ndim=1
		coeffs=reshape(coeffs,(1,length(coeffs)))
	end
	return ForceFromScipySplines(knots,coeffs,k,der,ndim)
end

function forceUpdate!(force::ForceFromScipySplines,f::Vector{TF}, x::Vector{TF})  where{TF<:AbstractFloat}
	for d in 1:force.ndim
		f[d]=scipy_interpolate.splev(x[1], (force.knots, force.coeffs_splines[d,:], force.k),force.der)[]
	end
	# for fix in force.fixes
	# 	apply_fix!(fix,x,f)
	# end
	return 0
end

function addFix(force::FP,fix::AF) where {FP<:AbstractForce,AF<:AbstractFix}
	#Ici on doit initil
	push!(force.fixes,fix)
end


#La fonction qui calcule la force a actualiser pour n'avoir que forceUpdate! from basis
function forceUpdate!(force::ForceFromPotential,f::Vector{TF}, x::Vector{TF})  where{TF<:AbstractFloat}
	force.gradV!(f,x)
	f.*=-1
	# for fix in force.fixes
	# 	apply_fix!(fix,x,f)
	# end
	return 0
end



function forceUpdate!(force::FB,f::Vector{TF}, x::Vector{TF})  where{FB<:AbstractForceFromBasis,TF<:AbstractFloat}
	for d in 1:force.ndim
		f[d]=force.basis[d](x[1]) # TODO, dealing with nd function
	end
	# for fix in force.fixes
	# 	apply_fix!(fix,x,f)
	# end
	return 0
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
