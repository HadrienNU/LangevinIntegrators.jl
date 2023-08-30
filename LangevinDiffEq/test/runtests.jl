using SafeTestsets
using Pkg

function activate_gpu_env()
    Pkg.activate("gpu")
    Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
    Pkg.instantiate()
end

const LONGER_TESTS = false

const GROUP = get(ENV, "GROUP", "All")

const is_APPVEYOR = Sys.iswindows() && haskey(ENV,"APPVEYOR")

  if !is_APPVEYOR && (GROUP == "All" || GROUP == "AlgConvergence")
    @time @safetestset "Dynamical SDE Tests" begin include("inertial.jl") end
  end
