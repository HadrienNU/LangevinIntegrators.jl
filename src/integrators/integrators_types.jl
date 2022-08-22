
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

function InitState(integrator::OverdampedIntegrator;dim=1)
    return InitState!(zeros(dim), integrator)
end

function InitState(integrator::OverdampedIntegrator, init_cond::AbstractInitCond ;  id=1)
    state = InitState(integrator)
    state.x=generate_initcond(init_cond; id=id)
    return state
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

function InitState(integrator::InertialIntegrator, init_cond_pos ::AbstractInitCond, init_cond_vitesse::AbstractInitCond ; id=1)
    state = InitState(integrator)
    state.x=generate_initcond(init_cond_pos; id=id)
    state.v=generate_initcond(init_cond_vitesse; id=id)
    return state
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

function InitState(integrator::InertialIntegrator, init_cond_pos ::AbstractInitCond, init_cond_vitesse::AbstractInitCond, init_cond_hidden::AbstractInitCond ;  id=1)
    state = InitState(integrator)
    state.x=generate_initcond(init_cond_pos; id=id)
    state.v=generate_initcond(init_cond_vitesse; id=id)
    state.h=generate_initcond(init_cond_hidden; dim=integrator.dim_tot-dim,id=id)
    return state
end
