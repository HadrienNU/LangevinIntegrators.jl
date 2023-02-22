module Plumed
    import LangevinIntegrators: AbstractFix, init_fix, apply_fix!, close_fix

    struct plumed{TF<: AbstractFloat} <: AbstractFix
            function plumed(plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1, delta_t::Float64=1e-3; temperature=1.0)
                return new()
            end
    end

    function init_fix(fix::plumed; kwargs...) # Différent du constructeur au sens que c'est appelé au début de chaque traj
        @warn("Plumed has not been loaded, this is an empty fix.")
        return fix
    end

    function apply_fix!(fix::plumed,x::Array{TF},f::Array{TF}; kwargs...) where {TF<:AbstractFloat}
    	return false
    end


    function close_fix(fix::plumed)
    end

    function addPlumed!(integrator::AbstractIntegrator,plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1; temperature=1.0)
        addFix!(integrator, plumed(plumed_input_file, plumed_log_file, dim, integrator.Δt; temperature=temperature))
    end

    export plumed
    export init_fix
    export apply_fix!
    export close_fix
    export addPlumed!

end
