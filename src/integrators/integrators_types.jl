
abstract type AbstractIntegrator end

abstract type OverdampedIntegrator <: AbstractIntegrator end
abstract type InertialIntegrator <: AbstractIntegrator end
abstract type HiddenIntegrator <: AbstractIntegrator end


abstract type AbstractState end

abstract type AbstractOverdampedState <: AbstractState end
abstract type AbstractOverdampedMemoryHiddenState <: AbstractOverdampedState end
abstract type AbstractInertialState <: AbstractState end
abstract type AbstractMemoryKernelState <: AbstractInertialState end
abstract type AbstractMemoryHiddenState <: AbstractInertialState end


# InitState from other State

function InitState!(s::AbstractOverdampedState, integrator::OverdampedIntegrator)
    return InitState!(s.x, integrator)
end

function InitState(s::AbstractOverdampedState, integrator::OverdampedIntegrator)
    return InitState(deepcopy(s.x), integrator)
end


function InitState!(s::AbstractInertialState, integrator::InertialIntegrator)
    return InitState!(s.x,s.v,integrator)
end

function InitState(s::AbstractInertialState, integrator::InertialIntegrator)
    return InitState(deepcopy(s.x),deepcopy(s.v), integrator)
end

function InitState!(s::AbstractMemoryHiddenState, integrator::HiddenIntegrator)
    return InitState!(s.x,s.v,s.h,integrator)
end

function InitState(s::AbstractMemoryHiddenState, integrator::HiddenIntegrator)
    return InitState(deepcopy(s.x),deepcopy(s.v),deepcopy(s.h), integrator)
end
