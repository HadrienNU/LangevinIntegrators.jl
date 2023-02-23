# LangevinIntegratorsPlumedExt

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HadrienNU.github.io/LangevinIntegrators.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HadrienNU.github.io/LangevinIntegrators.jl/dev/)
[![Build Status](https://github.com/HadrienNU/LangevinIntegrators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HadrienNU/LangevinIntegrators.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl/branch/main/graph/badge.svg?token=vlYbCnFhac)](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl)
<!-- [![Coverage](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl) -->

# Plumed integration for LangevinIntegrators

  The integrator can be coupled to Plumed (https://www.plumed.org/). Plumed should be set up a priori in your system.
  The library should be accessible at (pre)compile time, to test it run ldd <path to>/lib/libplumed.so  where <path to> is the path used as argument in the configure of plumed. Usually /usr/local. If that does not work add plumed library path to LD_LIBRARY_PATH.
  In general, it should work out of the box when plumed is installed in standard location.
  In case we need to specify location of plumed include file, we can use environment variable PLUMED_INCLUDE_PATH to setup the path to PLUMED_INCLUDE_PATH/plumed/wrapper/Plumed.h

  It have been tested with Plumed v2.5 to 2.7.


  As plumed consider to observe a 3D system of N atoms, the coordinate of the system are feeded to plumed as (for a 5D systems)

    (x_1,x_2,x_3) (x_4,x_5,0)

  Then coordinate can be obtained in plumed via the POSITION keyword, using ATOM=1 to get the first three coordinates in x,y,z, and so forth.


  Note that Plumed add a significant overhead to the run time. So depending of your need that could be a good move to implemented the needed functionnality into a Fix (see fix.jl for examples).
