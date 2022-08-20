"""
LangevinIntegrators.jl

Library to generate trajectories from (Generalized) Langevin Equation
"""

module LangevinIntegrators

using ConfParser
using TimerOutputs # Somes timer for debug #Gérer aussi le module de @timeit_debug

#Package for the force evaluation
using ForwardDiff # For automatic differentiation of potential
using ApproxFun # Various basis function
using BSplineKit # Bsplines function

const timer = TimerOutput() # A global timer to time computation when debugging

function similar(x::TF) where{TF<:AbstractFloat} #Helper function for the 1D case
    return x
end


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
