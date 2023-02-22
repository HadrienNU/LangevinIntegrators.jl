module Plumed
    using CBinding

    import LangevinIntegrators: AbstractFix, init_fix, apply_fix!, close_fix
    import LangevinIntegrators: AbstractIntegrator, addFix!

    # try
        path_plumed = get(ENV,"PLUMED_INCLUDE_PATH","/usr/local/include")

        c`-llibplumed.so -I$(path_plumed)`
        c"#include <plumed/wrapper/Plumed.h>"

        mutable struct plumed{TF<: AbstractFloat} <: AbstractFix

            plumedmain

            plumed_input_file::String
            plumed_log_file::String

            delta_t::Float64
            temperature::Float64

            dim::Int64
            natoms::Int64
            position_plumed::Vector{TF}
            forces_plumed::Vector{TF}
            #Some Ref to variables to pass and get back from plumed
            ref_bias::Ref{Cdouble}
            ref_stop_condition::Ref{Cint}
            ref_needs_energy::Ref{Cint}

            box::Array{Float64} # I don't really need those, but there are stored to not loose the reference
            virial::Array{Float64}
            masses::Array{Float64}

            bias_energy::Float64

        end

        function plumed(plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1, delta_t::Float64=1e-3; temperature=1.0)
            #Check that plumed_input_file exist
            isfile(plumed_input_file) || @warn "Plumed input file does not exist"

            natoms=div(dim,3)+(dim%3!=0)

            # We create all the needed array
            box=zeros(3,3)
            box[1,1]=1.0 # TODO: A changer et voir par quoi?
            box[2,2]=1.0
            box[3,3]=1.0

            return plumed(nothing,plumed_input_file,plumed_log_file,delta_t,temperature,dim,natoms,zeros(3*natoms),zeros(3*natoms),Ref{Cdouble}(0.0),Ref{Cint}(0),Ref{Cint}(0),box,zeros(3,3),ones(natoms),0.0)
        end


        function plumed_init!(plumed_struct; kwargs...)


            isnothing(plumed_struct.plumedmain) || plumed_finalize(plumed_struct) # If not nothing, finalize to reinitialize

            ptr_plumed_file=pointer(plumed_struct.plumed_input_file)
            ptr_plumed_log=pointer(plumed_struct.plumed_log_file)

            plumed_struct.plumedmain=c"plumed_create"()

            ref_plumed_api_version=Ref{Cint}(0)
            c"plumed_cmd"(plumed_struct.plumedmain,"getApiVersion",ref_plumed_api_version) #Do Some check on plumed_api (en particulier)
            plumed_api_version=ref_plumed_api_version[]

            get(kwargs,:verbose,0) > 1 && println("Plumed API VERSION $plumed_api_version")

            if plumed_api_version < 6 || plumed_api_version > 11  # TODO: Faire des Tests sur l'api de 6 à 9
                println("Incompatible API version $plumed_api_version for PLUMED. Only Plumed 2.5.x, 2.6.x, and 2.7.x are tested and supported.")
            end

            c"plumed_cmd"(plumed_struct.plumedmain,"setRealPrecision",Ref(8));     # Pass a pointer to an integer containing the size of a real number (4 or 8)
            c"plumed_cmd"(plumed_struct.plumedmain,"setMDEnergyUnits",Ref(1.0));        # Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
            c"plumed_cmd"(plumed_struct.plumedmain,"setMDLengthUnits",Ref(1.0));        # Pass a pointer to the conversion factor between the length unit used in your code and nm
            c"plumed_cmd"(plumed_struct.plumedmain,"setMDTimeUnits",Ref(1.0));            # Pass a pointer to the conversion factor between the time unit used in your code and ps

            # This is valid only if API VERSION > 3
            c"plumed_cmd"(plumed_struct.plumedmain,"setMDChargeUnits",Ref(1.0));        # Pass a pointer to the conversion factor between the charge unit used in your code and e
            c"plumed_cmd"(plumed_struct.plumedmain,"setMDMassUnits",Ref(1.0));            # Pass a pointer to the conversion factor between the mass unit used in your code and amu

            c"plumed_cmd"(plumed_struct.plumedmain,"setPlumedDat",ptr_plumed_file);            # Pass the name of the plumed input file from the md code to plumed
            c"plumed_cmd"(plumed_struct.plumedmain,"setLogFile",ptr_plumed_log);                     # Pass the file  on which to write out the plumed log (to be created)


            name_integrators="LangevinIntegrators"
            c"plumed_cmd"(plumed_struct.plumedmain,"setMDEngine",pointer(name_integrators));                # Pass the name of your md engine to plumed (now it is just a label)

            c"plumed_cmd"(plumed_struct.plumedmain,"setNatoms",Ref{Cint}(plumed_struct.natoms));                    # Pass a pointer to the number of atoms in the system to plumed
            c"plumed_cmd"(plumed_struct.plumedmain,"setTimestep",Ref(plumed_struct.delta_t));                 # Pass a pointer to the molecular dynamics timestep to plumed

            # This is valid only if API VERSION > 1
            c"plumed_cmd"(plumed_struct.plumedmain,"setKbT",Ref(1.0*plumed_struct.temperature));                          # Pointer to a real containing the value of kbT

            # This is valid only if API VERSION > 2
            c"plumed_cmd"(plumed_struct.plumedmain,"setRestart",Ref(0));                      # Pointer to an integer saying if we are restarting (zero means no, one means yes)

            # Calls to do the actual initialization (all the above commands must appear before this call)
            c"plumed_cmd"(plumed_struct.plumedmain,"init",C_NULL);                            # Do all the initialization of plumed


            c"plumed_cmd"(plumed_struct.plumedmain,"setBox",Ref(plumed_struct.box,1));                  # Pass a pointer to the first element in the box share array to plumed

            return plumed_struct
        end

        function plumed_step!(plumed_struct::plumed,n,positions,forces,energy)

            plumed_struct.position_plumed[1:plumed_struct.dim]=positions
            plumed_struct.forces_plumed[1:plumed_struct.dim]=forces
            if (plumed_struct.dim+1)<=3*plumed_struct.natoms
                plumed_struct.position_plumed[(plumed_struct.dim+1):3*plumed_struct.natoms].=0.0
                plumed_struct.forces_plumed[(plumed_struct.dim+1):3*plumed_struct.natoms].=0.0
            end

            c"plumed_cmd"(plumed_struct.plumedmain,"setStep",Ref{Cint}(n));                      # Pass a pointer to the current timestep to plumed
            c"plumed_cmd"(plumed_struct.plumedmain,"setStopFlag",plumed_struct.ref_stop_condition)
            c"plumed_cmd"(plumed_struct.plumedmain,"setMasses",Ref(plumed_struct.masses,1));                 # Pass a pointer to the first element in the masses array to plumed
            c"plumed_cmd"(plumed_struct.plumedmain,"setVirial",Ref(plumed_struct.virial,1));
            c"plumed_cmd"(plumed_struct.plumedmain,"setPositions",Ref(plumed_struct.position_plumed,1));
            c"plumed_cmd"(plumed_struct.plumedmain,"setForces",Ref(plumed_struct.forces_plumed,1));                 # Pass a pointer to the first element in the foces array to plumed


            c"plumed_cmd"(plumed_struct.plumedmain,"prepareCalc",C_NULL);

            c"plumed_cmd"(plumed_struct.plumedmain,"isEnergyNeeded",plumed_struct.ref_needs_energy);# assuming flag is an int, that will be set to 0 if energy is not needed and 1 if it is needed
            if plumed_struct.ref_needs_energy[] == 1
                c"plumed_cmd"(plumed_struct.plumedmain,"setEnergy",Ref{Cfloat}(energy));                  # Pass a pointer to the current value of the potential energy to plumed?
            end
            c"plumed_cmd"(plumed_struct.plumedmain,"performCalc",C_NULL);

            c"plumed_cmd"(plumed_struct.plumedmain,"getBias",plumed_struct.ref_bias)

            forces[1:plumed_struct.dim]=plumed_struct.forces_plumed[1:plumed_struct.dim]
            return plumed_struct.ref_stop_condition[], plumed_struct.ref_bias[]
        end


        function plumed_finalize(plumed_struct::plumed)
            c"plumed_finalize"(plumed_struct.plumedmain)
            plumed_struct.plumedmain=nothing
        end


        function init_fix(fix::plumed; kwargs...) # Différent du constructeur au sens que c'est appelé au début de chaque traj
            plumed_init!(fix)
        end

        function apply_fix!(fix::plumed,x::Array{TF},f::Array{TF}; kwargs...) where {TF<:AbstractFloat}
            step=get(kwargs,:step,1)
            stop_cond, energy = plumed_step!(fix,step,x,f,0.0) #We just compute the extra energy that the fix is bringing
            fix.bias_energy=energy
        	return stop_cond
        end


        function close_fix(fix::plumed)
            plumed_finalize(fix)
        end


    # catch e  # Allow to proceed even when plumed is not found
    #     @warn("Unable to load plumed extension: ",e)
    #     @warn("A empty plumed module have been loaded instead, it has not effect.")
    #
    #     struct plumed{TF<: AbstractFloat} <: AbstractFix
    #             function plumed(plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1, delta_t::Float64=1e-3; temperature=1.0)
    #                 return new()
    #             end
    #     end
    #
    #     function init_fix(fix::plumed; kwargs...) # Différent du constructeur au sens que c'est appelé au début de chaque traj
    #         @warn("Plumed has not been loaded, this is an empty fix.")
    #         return fix
    #     end
    #
    #     function apply_fix!(fix::plumed,x::Array{TF},f::Array{TF}; kwargs...) where {TF<:AbstractFloat}
    #     	return false
    #     end
    #
    #
    #     function close_fix(fix::plumed)
    #     end
    #
    #
    # end

    function addPlumed!(integrator::AbstractIntegrator,plumed_input_file::String, plumed_log_file::String="p.log", dim::Int64=1; temperature=1.0)
        addFix!(integrator, plumed(plumed_input_file, plumed_log_file, dim, integrator.Δt; temperature=temperature))
    end



    export plumed
    export init_fix
    export apply_fix!
    export close_fix
    export addPlumed!

end
