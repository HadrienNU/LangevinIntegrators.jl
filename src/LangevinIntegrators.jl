"""
LangevinIntegrators.jl

Library to generate trajectories from (Generalized) Langevin Equation
"""

module LangevinIntegrators

using ConfParser
# using TimerOutputs # Somes timer for debug #Gérer aussi le module de @timeit_debug
# En vrai il faudrait le remplacer par PkgBenchmark.jl ou juste @profile
# PkgBenchmark.jl ça marche avec BenchmarkTools.jl

#Package for the force evaluation
using ForwardDiff # For automatic differentiation of potential
using ApproxFun # Various basis function
using BSplineKit # Bsplines function

using LinearAlgebra

# const timer = TimerOutput() # A global timer to time computation when debugging




include("potentials.jl")
include("force.jl")
export ForceFromPotential
export ForceFromBasis

#Plumed module
include("plumed.jl")


#Integrators
include("integrators/integrators_types.jl")

include("integrators/overdamped_em.jl")
export EM
#include("integrators/bbk.jl")
#export BBK
#include("integrators/g_jf.jl")
#export GJF
#include("integrators/aboba.jl")
#export ABOBA
#include("integrators/baoab.jl")
#export BAOAB
#include("integrators/verlet.jl")
#export Verlet


include("params.jl")

export read_conf,LangevinParams

include("run.jl")
export init_trajectory
export run_trajectory!

include("initialization.jl")


end
