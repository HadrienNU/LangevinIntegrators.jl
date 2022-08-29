#=
Les functions pour l'analyse de la fin de la trajectoire. Donc plutôt du calcul de FPT
En vrai je pense qu'on va virer cette partie, il y a assez peu de cas d'usage, et on peut en faire l'essentiel via plumed.
Si la performance est une issue, alors on sauvegarde la trajectoire dans un fichier et on l'analyse avec plumed plus tard
=#

abstract type AbstractObserver end

abstract type AbstractStatisticalObs <: AbstractObserver end

"""
Function to initialize the observers
"""

mutable struct FPT_calc <: AbstractStatisticalObs
    reached::Array{Bool}
    fpt::Array{Float64} #Final time of the trajectories
    FPT_calc() = new(empty([], Bool), empty([], Float64))
end

function obs_end_traj!(observer::FPT_calc, state::AbstractState; end_time, stoppingCriterion)
    push!(observer.reached,(stoppingCriterion != 0))
    push!(observer.fpt,end_time)
end


function calc_MFPT(observer::FPT_calc)
    return sum(observer.fpt) / observer.reached
end

# function calcFPT_cdf(observer::FPT_calc) # Do some Kaplan Meier analysis, en vrai on juste utiliser le package Survival
#     body
# end
