#  On va metrre ici ce qu'il faut au fur et à mesure on étendra après


# Pour définir les équations, on pourra toujours faire une surcouche qui prends tous les parmètres en entrées comme actuellement et retourne ce qu'il faut

  # Es-ce qu'on leut donne des paramètres, genre gamma const?
struct BAOAB <: StochasticDiffEqAlgorithm
end

# Par exemple, mais c'est probablement de l'optimisation prémature
struct BAOAB{const} <: StochasticDiffEqAlgorithm end
BAOAB(const=true) = BAOAB{const}()  # Et comme ça on peut faire du dispatch sur les propriétés de l'intégrateur

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
GJF()=GJ("II")
