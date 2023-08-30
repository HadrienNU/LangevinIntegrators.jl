module LangevinDiffEq

using Reexport
@reexport using DiffEqBase
@reexport using StochasticDiffEq
using RecursiveArrayTools

import StochasticDiffEq: StochasticDiffEqAlgorithm, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,is_diagonal_noise

# Functions to be extanded
import StochasticDiffEq: alg_cache, initialize!, perform_step!, alg_compatible,@cache

using DiffEqBase: @..

using UnPack, MuladdMacro
# On pourrait aussi faire @importall(StochasticDiffEq), en utilisant la package ImportAll.jl

using Memoization

include("algorithms.jl")

export GJ, VEC
export ABOBA,BAOAB, OBABO


macro cache(expr)
  name = expr.args[2].args[1].args[1]
  fields = expr.args[3].args[2:2:end]
  cache_vars = Expr[]
  rand_vars = Expr[]
  jac_vars = Pair{Symbol,Expr}[]
  ratenoise_vars = Expr[]
  for x in fields
    if x.args[2] == :uType || x.args[2] == :rateType ||
       x.args[2] == :kType || x.args[2] == :uNoUnitsType #|| x.args[2] == :possibleRateType
      push!(cache_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :JCType
      push!(cache_vars,:(c.$(x.args[1]).duals...))
    elseif x.args[2] == :GCType
      push!(cache_vars,:(c.$(x.args[1]).duals))
    elseif x.args[2] == :DiffCacheType
      push!(cache_vars,:(c.$(x.args[1]).du))
      push!(cache_vars,:(c.$(x.args[1]).dual_du))
    elseif x.args[2] == :JType || x.args[2] == :WType
      push!(jac_vars,x.args[1] => :(c.$(x.args[1])))
    elseif x.args[2] == :randType
      push!(rand_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :rateNoiseType || x.args[2] == :rateNoiseCollectionType
      # Should be a pair for handling non-diagonal
      push!(ratenoise_vars,:(c.$(x.args[1])))
    end
  end
  quote
    $expr  # Faire une pull request pour avoir $(esc(expr))
    $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
    $(esc(:jac_iter))($(esc(:c))::$name) = tuple($(jac_vars...))
    $(esc(:rand_cache))($(esc(:c))::$name) = tuple($(rand_vars...))
    $(esc(:ratenoise_cache))($(esc(:c))::$name) = tuple($(ratenoise_vars...))
  end
end

include("caches/inertial_caches.jl")
include("caches/splitting_caches.jl")

include("perform_step/dynamical.jl")
include("perform_step/gj.jl")
include("perform_step/splitting.jl")


end # module LangevinDiffEq
