
  # Es-ce qu'on leut donne des paramètres, genre gamma const?
struct BAOAB <: StochasticDiffEqAlgorithm
end

# # Par exemple, mais c'est probablement de l'optimisation prémature
# struct BAOAB{const} <: StochasticDiffEqAlgorithm end
# BAOAB(const=true) = BAOAB{const}()  # Et comme ça on peut faire du dispatch sur les propriétés de l'intégrateur

struct ABOBA <: StochasticDiffEqAlgorithm
end

struct OBABO <: StochasticDiffEqAlgorithm
end

struct BBK <: StochasticDiffEqAlgorithm
end

struct VEC <: StochasticDiffEqAlgorithm
end

struct GJ <: StochasticDiffEqAlgorithm
  type::String
end
GJ(;type="I") = GJ(type)
# GJF()=GJ("II")
# Ou aors via GJ{I},etc..


alg_compatible(prob::DiffEqBase.AbstractSDEProblem,alg::GJ) = true  # Ca devrait être is split problem
