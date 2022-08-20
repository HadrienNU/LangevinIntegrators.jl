
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
