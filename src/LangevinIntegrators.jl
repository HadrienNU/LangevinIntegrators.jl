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
    @warn("Unable to load plumed extension: ",e)
    @warn("A empty plumed module have been loaded instead, it has not effect.")

    include("plumed_replacement.jl")

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
