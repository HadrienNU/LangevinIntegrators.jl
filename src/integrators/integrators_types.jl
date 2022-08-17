
abstract type AbstractIntegrator end

abstract type OverdampedIntegrator <: AbstractIntegrator end
abstract type InertialIntegrator <: AbstractIntegrator end


abstract type AbstractState end

abstract type AbstractOverdampedState <: AbstractState end
abstract type AbstractInertialState <: AbstractState end
abstract type AbstractMemoryKernelState <: AbstractState end
abstract type AbstractMemoryHiddenState <: AbstractState end

#Pas sur que ça soit necessaire, puisque je peux toujours récupérer les infos depuis l'état de l'intégrateur
#mutable struct OverdampedState{TF<:AbstractFloat} <:AbstractOverdampedState
#    x::Vector{TF}
#end

#mutable struct InertialState{TF<:AbstractFloat} <:AbstractInertialState
#    x::Vector{TF}
#	v::Vector{TF}
#end

