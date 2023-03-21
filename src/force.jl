
#S'il y a un object plumed il doit être stocké dans la force?
# En tout cas, force update doit appeler plumed_step
# La force doit stocker une liste de Fix, dont celui de plumed
# Function a créer

abstract type AbstractForce end

abstract type AbstractForceFromBasis <: AbstractForce end

function addFix!(force::FP, fix::AbstractFix) where {FP<:AbstractForce}
    push!(force.fixes, fix)
    return force
end

function addFix!(integrator::AbstractIntegrator, fix::AbstractFix) # Add it to the integrator
    addFix!(integrator.force, fix)
end

# TODO: add a function to get current potential energy
"""
ForceFromPotential(potential)

Set up force using the derivative of a potential.

### Fields

* potential   - Potential name, should be the name of a function implemented in LangevinIntegrators
"""
struct ForceFromPotential{G<:Function} <: AbstractForce
    cfg::ForwardDiff.GradientConfig
    V::G # On stocke aussi la valeur du potential
    ndim::Int
    fixes::Array{AbstractFix} # List of fix to apply
end


function ForceFromPotential(potential::Union{String,Function}, ndim = 1::Int; x₀::Union{Vector{TF},Nothing} = nothing, potential_kwargs::Dict = Dict()) where {TF<:AbstractFloat}
    if isnothing(x₀)
        x₀ = zeros(ndim) #Vector{Float64}(undef, ndim)
    end
    if !isempty(methods(potential))
        V = x -> potential(x, potential_kwargs...)
    else
        V = x -> getfield(LangevinIntegrators, Symbol(potential))(x,potential_kwargs...) # Eventuellement créer un sous module avec juste les potentials pour ne pas chercher partout
    end
    cfg = ForwardDiff.GradientConfig(V, x₀)
    return ForceFromPotential(cfg, V, ndim, Vector{AbstractFix}(undef, 0))
end


function forceUpdate!(
    force::ForceFromPotential,
    f::Vector{TF},
    x::Vector{TF};
    kwargs...,
) where {TF<:AbstractFloat}
    ForwardDiff.gradient!(f, force.V, x, force.cfg)
    f .*= -1
    stop_condition = false
    if get(kwargs, :applyFix, true)
        for fix in force.fixes
            out_stop_condition = apply_fix!(fix, x, f; kwargs...)
            stop_condition = stop_condition || out_stop_condition
        end
    end
    return stop_condition
end

"""
ForceFromBasis(potential)

Set up force using a set of basis function from ApproxFun package.
See https://juliaapproximation.github.io/ApproxFun.jl/latest/usage/spaces/ for a list of available basis

### Fields

* type   - Basis name, should be the name of a space implemented in ApproxFun
* coeffs - Coefficient of the force in the basis
"""
struct ForceFromBasis <: AbstractForceFromBasis #
    basis::Vector{Fun}
    ndim::Int
    fixes::Array{AbstractFix} # List of fix to apply
end

function ForceFromBasis(type::String, coeffs::Array{TF}) where {TF<:AbstractFloat}
    if ndims(coeffs) >= 2
        ndim = size(coeffs)[1]
    else
        ndim = 1
    end
    function_space = getfield(ApproxFun, Symbol(type))
    basis = Vector{Fun}(undef, ndim)
    for d = 1:ndim
        basis[d] = Fun(function_space, coeffs[d, :])
    end
    return ForceFromBasis(basis, ndim, Vector{AbstractFix}(undef, 0))
end

function forceUpdate!(
    force::FB,
    f::Vector{TF},
    x::Vector{TF},
    ;
    kwargs...,
) where {FB<:AbstractForceFromBasis,TF<:AbstractFloat}
    for d = 1:force.ndim
        f[d] = force.basis[d](x[1]) # TODO, dealing with nd function
    end
    stop_condition = false
    if get(kwargs, :applyFix, true)
        for fix in force.fixes
            out_stop_condition = apply_fix!(fix, x, f; kwargs...)
            stop_condition = stop_condition || out_stop_condition
        end
    end
    return stop_condition
end


"""
ForceFromSplines(k, knots, coeffs)

Set up force using splines from BSplineKit.

### Fields

* k   - Degree of the splines
* knots     - Knots
* coeffs     - Coefficients
"""
struct ForceFromSplines <: AbstractForceFromBasis # use BSplineKit # The struct is writed to have the same interface than Force from basis to have only one forceUpdate
    basis::Array#{Splines}
    ndim::Int
    fixes::Array{AbstractFix} # List of fix to apply
end

function ForceFromSplines(
    k::Int,
    knots::Array{TF},
    coeffs::Array{TF};
    der = 0,
) where {TF<:AbstractFloat}
    if ndims(coeffs) >= 2 #TODO change for a Vector of vector
        ndim = size(coeffs)[1]
        nb_coeffs = size(coeffs)[2]
    else
        ndim = 1
        nb_coeffs = length(coeffs)
        coeffs = reshape(coeffs, (1, nb_coeffs))
    end
    #Note the order of the splines is k+1 with k the degree
    B = BSplineBasis(BSplineOrder(k + 1), knots; augment = Val(false))
    basis = Vector{SplineExtrapolation}(undef, ndim)
    for d = 1:ndim
        basis[d] = BSplineKit.SplineExtrapolations.extrapolate(
            BSplineKit.Derivative(der) * Spline(
                B,
                coeffs[d, 1:(nb_coeffs == length(knots) ? nb_coeffs - (k + 1) : nb_coeffs)],
            ),
            BSplineKit.SplineExtrapolations.Smooth(),
        )
    end
    return ForceFromSplines(basis, ndim, Vector{AbstractFix}(undef, 0))
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
"""
ForceFromScipySplines(k, knots, coeffs)

Set up force using splines from scipy.
This is slower than julia implemenation of splines but allow to use the python implementation if needed.

### Fields

* k   - Degree of the splines
* knots     - Knots
* coeffs     - Coefficients
"""
struct ForceFromScipySplines{TF<:AbstractFloat} <: AbstractForceFromBasis  # use Scipy
    knots::Array{TF}
    coeffs_splines::Array{TF}
    k::Int
    der::Int
    ndim::Int
    fixes::Array{AbstractFix} # List of fix to apply
end
function ForceFromScipySplines(
    k::Int,
    knots::Array{TF},
    coeffs::Array{TF};
    der = 0,
) where {TF<:AbstractFloat}
    if ndims(coeffs) >= 2 #TODO change for a Vector of vector
        ndim = size(coeffs)[1]
    else
        ndim = 1
        coeffs = reshape(coeffs, (1, length(coeffs)))
    end
    return ForceFromScipySplines(knots, coeffs, k, der, ndim, Vector{AbstractFix}(undef, 0))
end

function forceUpdate!(
    force::ForceFromScipySplines,
    f::Vector{TF},
    x::Vector{TF},
    ;
    kwargs...,
) where {TF<:AbstractFloat}
    for d = 1:force.ndim
        pylock() do # Lock Thread
            f[d] = scipy_interpolate.splev(
                x[1],
                (force.knots, force.coeffs_splines[d, :], force.k),
                force.der,
            )[]
        end
    end
    stop_condition = false
    if get(kwargs, :applyFix, true)
        for fix in force.fixes
            out_stop_condition = apply_fix!(fix, x, f; kwargs...)
            stop_condition = stop_condition || out_stop_condition
        end
    end
    return stop_condition
end



function forceUpdate(
    force::FP,
    x::Vector{TF};
    applyFix = false,
) where {FP<:AbstractForce,TF<:AbstractFloat}
    f = similar(x)
    forceUpdate!(force, f, x; applyFix = applyFix) # For the initialization or the computation of the force, do not include fix
    return f
end

function force_eval(force::FP, x::Vector{TF}) where {FP<:AbstractForce,TF<:AbstractFloat} # When we just want a single point evaluation
    return force_eval(force, [x])[1]
end

# A function to plot value of the force for comparison with python code
function force_eval(
    force::FP,
    x::Vector{Vector{TF}};
    applyFix = true,
) where {FP<:AbstractForce,TF<:AbstractFloat} # Evaluate over a bunch of point
    f = similar(x)
    for (n, val) in enumerate(x)
        f[n] = forceUpdate(force, val; applyFix = applyFix)
    end
    return f
end
