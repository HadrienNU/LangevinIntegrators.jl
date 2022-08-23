
abstract type AbstractIntegrator end

abstract type OverdampedIntegrator <: AbstractIntegrator end
abstract type InertialIntegrator <: AbstractIntegrator end
abstract type HiddenIntegrator <: AbstractIntegrator end
abstract type KernelIntegrator <: AbstractIntegrator end


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

function InitState(integrator::OverdampedIntegrator;dim=1)
    return InitState!(zeros(dim), integrator)
end

function InitState!(s::AbstractInertialState, integrator::InertialIntegrator)
    return InitState!(s.x,s.v,integrator)
end

function InitState(s::AbstractInertialState, integrator::InertialIntegrator)
    return InitState(deepcopy(s.x),deepcopy(s.v), integrator)
end

function InitState(integrator::InertialIntegrator;dim=1)
    return InitState!(zeros(dim),zeros(dim), integrator)
end

function InitState!(s::AbstractMemoryHiddenState, integrator::HiddenIntegrator)
    return InitState!(s.x,s.v,s.h,integrator)
end

function InitState(s::AbstractMemoryHiddenState, integrator::HiddenIntegrator)
    return InitState(deepcopy(s.x),deepcopy(s.v),deepcopy(s.h), integrator)
end

function InitState(integrator::HiddenIntegrator;dim=1)
    return InitState!(zeros(dim),zeros(dim),zeros(integrator.dim_tot-dim), integrator)
end

#InitState from init_cond generator

function InitState(integrator::OverdampedIntegrator, init_cond::Array ;  id=1)
    return InitState(generate_initcond(init_cond[1]; id=id),integrator)
end

function InitState(integrator::InertialIntegrator, init_cond::Array ; id=1)
    return InitState(generate_initcond(init_cond[1]; id=id),generate_initcond(init_cond[2]; id=id),integrator)
end

function InitState(integrator::HiddenIntegrator,init_cond::Array;  id=1)
    return InitState(generate_initcond(init_cond[1]; id=id),generate_initcond(init_cond[2]; id=id),generate_initcond(init_cond[3]; id=id),integrator)
end
