
  # Es-ce qu'on leut donne des paramètres, genre gamma const?
struct BAOAB{T} <: StochasticDiffEqAlgorithm
    kT::T
end
BAOAB(kT=1.0) = BAOAB(kT)

# # Par exemple, mais c'est probablement de l'optimisation prémature
# struct BAOAB{const} <: StochasticDiffEqAlgorithm end
# BAOAB(const=true) = BAOAB{const}()  # Et comme ça on peut faire du dispatch sur les propriétés de l'intégrateur

struct ABOBA{T} <: StochasticDiffEqAlgorithm
    kT::T
end
ABOBA(kT=1.0) = ABOBA(kT)

struct OBABO{T} <: StochasticDiffEqAlgorithm
    kT::T
end
OBABO(kT=1.0) = OBABO(kT)

struct BBK{T} <: StochasticDiffEqAlgorithm
    kT::T
end
BBK(kT=1.0) = BBK(kT)

struct VEC{T} <: StochasticDiffEqAlgorithm
    kT::T
end
VEC(kT=1.0) = VEC(kT)

struct GJ{T} <: StochasticDiffEqAlgorithm
  kT::T
  type::String
end
GJ(kT=1.0; type="I") = GJ(kT,type)
# GJF()=GJ("II")
# Ou aors via GJ{I},etc..


alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::GJ) = is_diagonal_noise(prob) # To be relased later
alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::VEC) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::BAOAB) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::OBABO) = is_diagonal_noise(prob)
alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::ABOBA) = is_diagonal_noise(prob)
