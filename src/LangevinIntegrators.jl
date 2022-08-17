"""
LangevinIntegrators.jl

Library to generate trajectories from (Generalized) Langevin Equation
"""

module LangevinIntegrators

using ConfParser
using ForwardDiff # For automatic differentiation of potential
using TimerOutputs
using ApproxFun

#Gérer aussi le module de @timeit_debug

include("potentials.jl")
include("force.jl")
export ForceFromPotential
export ForceFromBasis
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
