using CBinding

#TODO: ça a l'air compliqu" de faire ça en dehors du scope global, à voir
#Modifier les dossiers d'inclusion selon les variables d'environment
#En vrai, il faut surtout vérifier si les fichiers exists
pathlib="/usr/local/lib/"
pathinclude="/usr/local/include"

c`-L$(pathlib) -llibplumed.so -I$(pathinclude)`
c"""
  #include <plumed/wrapper/Plumed.h>
"""

struct Plumed{TF<: AbstractFloat}
    plumedmain
    dim::Int64
    natoms::Int64
    position_plumed::Vector{TF}
    forces_plumed::Vector{TF}
    #Some Ref to variables to pass and get back from plumed
    ref_energy::Ref{Cfloat}
    ref_stop_condition::Ref{Cint}
    ref_needs_energy::Ref{Cint}

    box::Array{Float64} # I don't really need those, but there are stored to not loose the reference
    virial::Array{Float64}
    masses::Array{Float64}
end

function Plumed(plumed_input_file,plumed_log_file="p.log",dim=1,delta_t=1e-3,temperture=1.0;path_to_plumed="/usr/local/")


    #Je dois récupérer une structure de param avec dim, delta_t, temperature, masses
    ptr_plumed_file=pointer(plumed_input_file)
    ptr_plumed_log=pointer(plumed_log_file)

    natoms=div(dim,3)+(dim%3!=0)

    # We create all the needed array
    box=zeros(3,3)
    box[1,1]=1.0 # TODO: A changer et voir par quoi?
    box[2,2]=1.0
    box[3,3]=1.0

    plumedmain=c"plumed_create"()

    ref_plumed_api_version=Ref{Cint}(0)
    c"plumed_cmd"(plumedmain,"getApiVersion",ref_plumed_api_version) #Do Some check on plumed_api (en particulier)
    plumed_api_version=ref_plumed_api_version[]
    println("API VERSION $plumed_api_version")
    if plumed_api_version < 5 || plumed_api_version > 10  # TODO: Faire des Tests sur l'api de 6 à 9
        println("Incompatible API version $plumed_api_version for PLUMED. Only Plumed 2.4.x, 2.5.x, and 2.6.x are tested and supported.")
    end

    c"plumed_cmd"(plumedmain,"setRealPrecision",Ref(8));     # Pass a pointer to an integer containing the size of a real number (4 or 8)
    c"plumed_cmd"(plumedmain,"setMDEnergyUnits",Ref(1.0));        # Pass a pointer to the conversion factor between the energy unit used in your code and kJ mol-1
    c"plumed_cmd"(plumedmain,"setMDLengthUnits",Ref(1.0));        # Pass a pointer to the conversion factor between the length unit used in your code and nm
    c"plumed_cmd"(plumedmain,"setMDTimeUnits",Ref(1.0));            # Pass a pointer to the conversion factor between the time unit used in your code and ps

    # This is valid only if API VERSION > 3
    c"plumed_cmd"(plumedmain,"setMDChargeUnits",Ref(1.0));        # Pass a pointer to the conversion factor between the charge unit used in your code and e
    c"plumed_cmd"(plumedmain,"setMDMassUnits",Ref(1.0));            # Pass a pointer to the conversion factor between the mass unit used in your code and amu

    c"plumed_cmd"(plumedmain,"setPlumedDat",ptr_plumed_file);            # Pass the name of the plumed input file from the md code to plumed
    c"plumed_cmd"(plumedmain,"setLogFile",ptr_plumed_log);                     # Pass the file  on which to write out the plumed log (to be created)


    name_integrators="LangevinIntegrators"
    c"plumed_cmd"(plumedmain,"setMDEngine",pointer(name_integrators));                # Pass the name of your md engine to plumed (now it is just a label)

    c"plumed_cmd"(plumedmain,"setNatoms",Ref{Cint}(natoms));                    # Pass a pointer to the number of atoms in the system to plumed
    c"plumed_cmd"(plumedmain,"setTimestep",Ref(delta_t));                 # Pass a pointer to the molecular dynamics timestep to plumed

    # This is valid only if API VERSION > 1
    c"plumed_cmd"(plumedmain,"setKbT",Ref(1.0*temperature));                          # Pointer to a real containing the value of kbT

    # This is valid only if API VERSION > 2
    c"plumed_cmd"(plumedmain,"setRestart",Ref(0));                      # Pointer to an integer saying if we are restarting (zero means no, one means yes)

    # Calls to do the actual initialization (all the above commands must appear before this call)
    c"plumed_cmd"(plumedmain,"init",C_NULL);                            # Do all the initialization of plumed


    c"plumed_cmd"(plumedmain,"setBox",Ref(box,1));                  # Pass a pointer to the first element in the box share array to plumed


    return Plumed(plumedmain,dim,natoms,zeros(3*natoms),zeros(3*natoms),Ref{Cfloat}(0.0),Ref{Cint}(0),Ref{Cint}(0),box,zeros(3,3),ones(natoms))
end


function plumed_step(plumed_struct::Plumed,n,positions,forces,energy)

    plumed_struct.position_plumed[1:plumed_struct.dim]=positions
    plumed_struct.forces_plumed[1:plumed_struct.dim]=forces
    if (plumed_struct.dim+1)<=3*plumed_struct.natoms
        plumed_struct.position_plumed[(dim+1):3*plumed_struct.natoms].=0.0
        plumed_struct.forces_plumed[(dim+1):3*plumed_struct.natoms].=0.0
    end

    c"plumed_cmd"(plumed_struct.plumedmain,"setStep",Ref{Cint}(n));                      # Pass a pointer to the current timestep to plumed
    c"plumed_cmd"(plumed_struct.plumedmain,"setStopFlag",plumed_struct.ref_stop_condition)
    c"plumed_cmd"(plumed_struct.plumedmain,"setMasses",Ref(plumed_struct.masses,1));                 # Pass a pointer to the first element in the masses array to plumed
    c"plumed_cmd"(plumed_struct.plumedmain,"setVirial",Ref(plumed_struct.virial,1));

    c"plumed_cmd"(plumed_struct.plumedmain,"setPositions",Ref(plumed_struct.position_plumed,1));
    c"plumed_cmd"(plumed_struct.plumedmain,"setForces",Ref(plumed_struct.forces_plumed,1));                 # Pass a pointer to the first element in the foces array to plumed


    c"plumed_cmd"(plumed_struct.plumedmain,"prepareCalc",C_NULL);

    c"plumed_cmd"(plumed_struct.plumedmain,"isEnergyNeeded",plumed_struct.ref_needs_energy);# assuming flag is an int, that will be set to 0 if energy is not needed and 1 if it is needed
    plumed_struct.ref_energy[]=energy
    if plumed_struct.ref_needs_energy[] == 1
        c"plumed_cmd"(plumed_struct.plumedmain,"setEnergy",plumed_struct.ref_energy);                  # Pass a pointer to the current value of the potential energy to plumed?
    end
    c"plumed_cmd"(plumed_struct.plumedmain,"performCalc",C_NULL);
    energy=plumed_struct.ref_energy[]

    println("Stop conditions ",plumed_struct.ref_stop_condition[], " Need energy ", plumed_struct.ref_needs_energy[])
    forces[1:plumed_struct.dim]=plumed_struct.forces_plumed[1:plumed_struct.dim]
    return plumed_struct.ref_stop_condition[], energy
end


function plumed_finalize(plumed_struct)
    c"plumed_finalize"(plumed_struct.plumedmain)
end

#Quelques infos à changer
delta_t=2e-3
dim=4
temperature=1.0

# A passer, en vrai il faudrait plutôt fournir pour l'energie une fonction callback qui ne calculera l'energie que si nécessaire.
energy=0.0
pos = [1.0,2.0,3.0,4.0] #Vector{Float}(undef,dim)
forces = [0.1,1.5,2.5,3.5] #Vector{Float}(undef,dim)

plumed_struct=Plumed("../examples/plumed.dat","../examples/plumed.log",dim,delta_t,temperature)

for n in 1:10
    plumed_step(plumed_struct,n,pos,forces,energy)
    println(forces)
end

plumed_finalize(plumed_struct)
