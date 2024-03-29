# LangevinIntegrators

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HadrienNU.github.io/LangevinIntegrators.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HadrienNU.github.io/LangevinIntegrators.jl/dev/)
[![Build Status](https://github.com/HadrienNU/LangevinIntegrators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HadrienNU/LangevinIntegrators.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl/branch/main/graph/badge.svg?token=vlYbCnFhac)](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl)
<!-- [![Coverage](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl) -->


A package to generate trajectories from a (generalized) Langevin equation. This a set of stochastic integrator in Julia of the Langevin equation.


# Overview


  This code is inspired by NDsimulator (https://ndsimulator.rtfd.io/) and the integrator are adapted from BasisMD.jl (https://github.com/gideonsimpson/BasicMD.jl)

# Setting up the integrator

  See the example folder for input file format.

    [space]
    ndim=1
    boundary=pbc
    low_limits=-10
    upp_limits=10

    [sampling]
    ntrajs=5
    dt=1e-3
    nsteps=5000
    integration=EM ; Integrateur to use, this is an overdamped integrator

    ; All parameters related to the setting of the force and the
    [physics]
    potential=Harmonic
    temperature=1.0

    ; How to generate initial conditions for the position
    [init_position]
    type=gaussian
    mean=0.0
    std=1.0

    [dump] ; How to output the trajectories
    n_save_iters=50  ; Frequency of save
    save_filename_pattern = trajectory_*.dat ; * are replaced by the id of the traj



  LangevinIntegrators is able to get parameters from the GLE_AnalysisEM package (https://github.com/HadrienNU/GLE_AnalysisEM)

# Structure of the code

  The flow is as follow, first define a force, via a potential or a basis. If wanted, you can add some fixes to the force to either communicate with plumed, add some extra forces or set tup a stop condition.
  Set up the integrator with this force and other physical parameters.

# Plumed integration

  The integrator can be coupled to Plumed (https://www.plumed.org/). Plumed should be set up a priori in your system.
  The library should be accessible at (pre)compile time, to test it run ldd <path to>/lib/libplumed.so  where <path to> is the path used as argument in the configure of plumed. Usually /usr/local. If that does not work add plumed library path to LD_LIBRARY_PATH.
  In general, it should work out of the box when plumed is installed in standard location.
  In case we need to specify location of plumed include file, we can use environment variable PLUMED_INCLUDE_PATH to setup the path to PLUMED_INCLUDE_PATH/plumed/wrapper/Plumed.h

  It have been tested with Plumed v2.5 to 2.7.


  As plumed consider to observe a 3D system of N atoms, the coordinate of the system are feeded to plumed as (for a 5D systems)

    (x_1,x_2,x_3) (x_4,x_5,0)

  Then coordinate can be obtained in plumed via the POSITION keyword, using ATOM=1 to get the first three coordinates in x,y,z, and so forth.


  Note that Plumed add a significant overhead to the run time. So depending of your need that could be a good move to implemented the needed functionnality into a Fix (see fix.jl for examples) or an observers (see observers.jl).


# Extend LangevinIntegrators.jl

  ## Extra integrators

  Integrators are in integrators folder. Two struct need to be set up, one for the integrator physical parameters (and that should include a force), that is more a immutable struct and one for current state of the system adpated to the integrator.

  Then evolution in time is provided by

    UpdateState!(state, integrator; kwargs...)

  Initialisation is done via InitState, that take a vector of AbstractInitCond

  ## Extra initial conditions

  Generator of initial conditions are in generate_init_cond.jl, they can be added to get_init_conditions function for automatic conversion between Dict value and struct

  ## Extra potentials

  Simply define a function that take the position as argument.

  ## Extra fix

  Fix derive from AbstractFix. There should be typically implemented when wanted function is not in plumed or you want to speed up code and then implment things that are in plumed such as the committor function.

    function init_fix(fix::AbstractFix; kwargs...) #Called at the start of each trajectory
    end

    function apply_fix!(fix::AbstractFix, x::Array{TF}, f::Array{TF}; kwargs...) where {TF<:AbstractFloat}
        return 0 # Should return a stop condition, 0 -> continue; if !=0 -> stop iteration
    end

    function close_fix(fix::AbstractFix) #Called at the end of each trajectory
    end

  Fixes are not called in the initialization of the force for the first timestep, except for ABOBA type integrator.


# TODO

  - Pass to run trajectories, un vector of struct of AbstractInitCond instead of the vector of Dict with params
  - Implement some extra space, such as a sphere
