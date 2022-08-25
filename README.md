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




  LangevinIntegrators is able to get parameters from the GLE_AnalysisEM package (https://github.com/HadrienNU/GLE_AnalysisEM)

# Plumed integration

  The integrator can be coupled to Plumed (https://www.plumed.org/). Plumed should be set up a priori in your system.
  The library should be accessible at (pre)compile time, to test it run ldd <path to>/lib/libplumed.so  where <path to> is the path used as argument in the configure of plumed. Usually /usr/local. If that does not work add plumed library path to LD_LIBRARY_PATH.
  In general, it should work out of the box when plumed is installed in standard location.
  In case we need to specify location of plumed include file, we can use environment variable PLUMED_INCLUDE_PATH to setup the path to PLUMED_INCLUDE_PATH/plumed/wrapper/Plumed.h

  It have been tested with Plumed v2.7 only, but that should work with version above 2.5.

  Note that Plumed add a significant overhead to the run time. So depending of your need that could be a good move to implemented the needed functionnality into a Fix (see fix.jl for examples) or an observers (see observers.jl).

# TODO
