using CBinding

#Modifier les dossiers d'inclusion selon les variables d'environment
c`-L/usr/local/lib/ -llibplumed.so -I/usr/local/include`

c"""
  #include <plumed/wrapper/Plumed.h>
"""
#Quelques infos à changer
plumed_input_file="plumed.dat"
plumed_log_file="plumed.log"
delta_t=2e-3
dim=4
natoms=div(dim,3)+(dim%3!=0)
temperature=1.0
pos = [1.0,2.0,3.0,4.0] #Vector{Float}(undef,dim)
forces = [0.1,1.5,2.5,3.5] #Vector{Float}(undef,dim)
position_plumed=Vector{Float64}(undef,3*natoms)
position_plumed[1:dim]=pos
position_plumed[(dim+1):3*natoms].=0.0
forces_plumed=zeros(3*natoms)
forces_plumed[1:dim]=forces
forces_plumed[(dim+1):3*natoms].=0.0
masses=ones(natoms)
virial=zeros(3,3)
box=zeros(3,3)
box[1,1]=1.0
box[2,2]=1.0
box[3,3]=1.0
#some Reference, I need to save the reference where I have want to access the data back
ref_api=Ref{Cint}(0) #Reference vers un Int  #Transformer en Ref{Cint}(0)
ref_natoms=Ref{Cint}(natoms)
ref_dt=Ref(delta_t)

ref_masses=Ref(masses,1)
ref_pos=Ref(position_plumed,1)
ref_forces=Ref(forces_plumed,1)
ref_virial=Ref(virial,1)
ref_box=Ref(box,1)

ref_stopcondition=Ref{Cbool}(false) #Reference vers un Int  #Transformer en Ref{Cint}(0)

ref_needs_energy=Ref{Cint}(0)

ptr_plumed_file=pointer(plumed_input_file)#Base.cconvert(Cstring,plumed_input_file)#Ref(plumed_input_file)
ptr_plumed_log=pointer(plumed_log_file)#Base.cconvert(Cstring,plumed_log_file)#Ref(plumed_log_file)

plumedmain=c"plumed_create"()
#Ensuite on a plus qu'à faire des appel comme ça, trouver comment passer des pointeurs

c"plumed_cmd"(plumedmain,"getApiVersion",ref_api) #Do Some check on plumed_api (en particulier)
plumed_api_version=ref_api[]
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
# c"plumed_cmd"(plumedmain,"setMPIComm",&MPI_COMM_WORLD);           # Pass a pointer to the MPI communicator to plumed
# notice that from fortran the command "setMPIFComm" should be used instead
c"plumed_cmd"(plumedmain,"setNatoms",ref_natoms);                    # Pass a pointer to the number of atoms in the system to plumed
c"plumed_cmd"(plumedmain,"setMDEngine",pointer("LangevinIntegrators"));                # Pass the name of your md engine to plumed (now it is just a label)
c"plumed_cmd"(plumedmain,"setLogFile",ptr_plumed_log);                     # Pass the file  on which to write out the plumed log (to be created)
c"plumed_cmd"(plumedmain,"setTimestep",ref_dt);                 # Pass a pointer to the molecular dynamics timestep to plumed

# This is valid only if API VERSION > 1
c"plumed_cmd"(plumedmain,"setKbT",Ref(1.0*temperature));                          # Pointer to a real containing the value of kbT

# This is valid only if API VERSION > 2
c"plumed_cmd"(plumedmain,"setRestart",Ref(0));                      # Pointer to an integer saying if we are restarting (zero means no, one means yes)

# Calls to do the actual initialization (all the above commands must appear before this call)
c"plumed_cmd"(plumedmain,"init",C_NULL);                            # Do all the initialization of plumed


#On ne les remplis qu'une seule fois pouisque ça ne change pas
c"plumed_cmd"(plumedmain,"setBox",ref_box);                  # Pass a pointer to the first element in the box share array to plumed

for n in 1:10

    forces_plumed[1:dim]=forces
    forces_plumed[(dim+1):3*natoms].=0.0
# Ensuite, c'est ce qui va dans la boucle
c"plumed_cmd"(plumedmain,"setStopFlag",ref_stopcondition)

c"plumed_cmd"(plumedmain,"setStep",Ref{Cint}(n));                      # Pass a pointer to the current timestep to plumed
c"plumed_cmd"(plumedmain,"setMasses",ref_masses);                 # Pass a pointer to the first element in the masses array to plumed
c"plumed_cmd"(plumedmain,"setPositions",ref_pos);
# c"plumed_cmd"(plumedmain,"setEnergy",&poteng);                  # Pass a pointer to the current value of the potential energy to plumed?
c"plumed_cmd"(plumedmain,"setForces",ref_forces);                 # Pass a pointer to the first element in the foces array to plumed
c"plumed_cmd"(plumedmain,"setVirial",ref_virial);

c"plumed_cmd"(plumedmain,"prepareCalc",C_NULL);

c"plumed_cmd"(plumedmain,"isEnergyNeeded",ref_needs_energy);# assuming flag is an int, that will be set to 0 if energy is not needed and 1 if it is needed
c"plumed_cmd"(plumedmain,"performCalc",C_NULL);
println("Stop conditions ",ref_stopcondition[], " Need energy ", ref_needs_energy[])
println(forces)
println(forces_plumed) # Les forces sont bien dans forces_plumed, il faut voir pour les retransférer dans forces
end

c"plumed_finalize"(plumedmain)
