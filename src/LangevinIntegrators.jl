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

using LinearAlgebra

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

include("plumed.jl") # Find a way to make this optionnal, using Require?



include("generate_initcond.jl")
export initialize_initcond

include("integrators/boundary_conditions.jl")
export SeparateSpace
export noBC, PBC, ReflectingBC


include("integrators/overdamped_em.jl")
export EM

#Inertial
include("integrators/bbk.jl")
export BBK
include("integrators/g_jf.jl")
export GJF
include("integrators/aboba.jl")
export ABOBA
include("integrators/baoab.jl")
export BAOAB
include("integrators/verlet.jl")
export Verlet

#Hidden
include("integrators/hidden_em.jl")
export EM_Hidden
include("integrators/hidden_aboba.jl")
export ABOBA_Hidden

#MemoryKernel


export UpdateState!

include("params.jl")

export set_from_conf, set_hidden_from_npz
export read_conf, read_integrator_conf, read_integrator_hidden_npz
export TrajsParams

include("save.jl")
export TrajectorySave, TrajectorySaveInertial, TrajectorySaveOnlyHidden, TrajectorySaveHidden


include("run.jl")
export generate_initial_conditions
export run_trajectory!
export run_trajectories
export run_fpt

end
