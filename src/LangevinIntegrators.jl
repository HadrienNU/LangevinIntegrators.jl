"""
LangevinIntegrators.jl

Library to generate trajectories from (Generalized) Langevin Equation
"""

module LangevinIntegrators

using ConfParser
using PyCall

const np = PyNULL()
const scipy_interpolate = PyNULL()

function __init__()
    copy!(np, pyimport("numpy"))
    copy!(scipy_interpolate, pyimport("scipy.interpolate"))
end

#Package for the force evaluation
using ForwardDiff # For automatic differentiation of potential
using ApproxFun # Various basis function
using BSplineKit # Bsplines function

using StaticArrays
using DataStructures

using LinearAlgebra
using FFTW

using DelimitedFiles

#Integrators
include("integrators/integrators_types.jl")
export InitState
export InitState!

include("fixes.jl")

export LWall, UWall, Quadratic

include("potentials.jl")
include("force.jl")
export ForceFromPotential
export ForceFromBasis
export ForceFromSplines
export ForceFromScipySplines

export force_eval

export addFix!

try
    include("plumed.jl") # Find a way to make this optionnal, using Require?
catch e  # Allow to proceed even when plumed is not found
    warn("Unable to load plumed extension: ",e)
    warn("A empty plumed module have been loaded instead, it has not effect.")

    module Plumed
        struct plumed{TF<: AbstractFloat} <: AbstractFix
                function plumed(plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1, delta_t::Float64=1e-3; temperature=1.0)
                    return new()
                end
        end
        function init_fix(fix::plumed; kwargs...) # Différent du constructeur au sens que c'est appelé au début de chaque traj
            warn("Plumed has not been loaded, this is an empty fix.")
            return fix
        end

        function apply_fix!(fix::plumed,x::Array{TF},f::Array{TF}; kwargs...) where {TF<:AbstractFloat}
        	return false
        end


        function close_fix(fix::plumed)
        end

        function addPlumed!(integrator::AbstractIntegrator,plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1; temperature=1.0)
            addFix!(integrator, plumed(plumed_input_file, plumed_log_file, dim, integrator.Δt; temperature=temperature))
        end
    end
end

include("generate_initcond.jl")
export initialize_initcond

include("integrators/boundary_conditions.jl")
export SeparateSpace
export noBC, PBC, ReflectingBC

include("integrators/correlated_noise.jl")


include("integrators/overdamped/overdamped_em.jl")
export EM

#Inertial
include("integrators/inertial/velocity_verlet.jl")
export VelocityVerlet
include("integrators/inertial/bbk.jl")
export BBK
include("integrators/inertial/splitting.jl")
export BAOAB
export OBABO
include("integrators/inertial/gj.jl")
export GJF
export GJ


include("integrators/inertial/position_verlet.jl")
export PositionVerlet
export ABOBA



#Hidden
include("integrators/hidden/hidden_em.jl")
export EM_Hidden
include("integrators/hidden/hidden_aboba.jl")
export ABOBA_Hidden

#MemoryKernel
include("integrators/kernel/kernel_em.jl")
export EM_Kernel
include("integrators/kernel/kernel_bbk.jl")
export BBK_Kernel
include("integrators/kernel/kernel_g_jf.jl")
export GJF_Kernel


export UpdateState!

include("params.jl")

export set_from_conf, set_hidden_from_npz
export read_conf, read_integrator_conf, read_integrator_hidden_npz
export TrajsParams

include("save.jl")
export TrajectorySave, TrajectorySaveInertial, TrajectorySaveOnlyHidden, TrajectorySaveHidden
export TransitionObserver

include("run.jl")
export generate_initial_conditions
export run_trajectory!
export run_trajectories
export run_fpt, run_transitions

end
