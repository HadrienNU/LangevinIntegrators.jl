"""
LangevinIntegrators.jl

Library to generate trajectories from (Generalized) Langevin Equation
"""

module LangevinIntegrators

using ConfParser
using NPZ
# En vrai il faudrait le remplacer par PkgBenchmark.jl ou juste @profile
# PkgBenchmark.jl ça marche avec BenchmarkTools.jl

#Package for the force evaluation
using ForwardDiff # For automatic differentiation of potential
using ApproxFun # Various basis function
using BSplineKit # Bsplines function

using LinearAlgebra

# const timer = TimerOutput() # A global timer to time computation when debugging


include("modifiers.jl")


include("potentials.jl")
include("force.jl")
export ForceFromPotential
export ForceFromBasis

export force_eval

#Plumed module
# include("plumed.jl")


include("generate_initcond.jl")

#Integrators
include("integrators/integrators_types.jl")
export InitState
export InitState!

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


include("observers.jl")

include("params.jl")

export read_conf,LangevinParams




include("run.jl")
export run_trajectory!
export run_trajectories_parallel


end
