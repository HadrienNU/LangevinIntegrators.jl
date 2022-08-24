using BenchmarkTools
using LangevinIntegrators
# using Random

#Defining somes benchmark
# J'ai besoin de tester, l'évaluation d'une trajectoire pour les différents intégrateurs

init_conds_args=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"Cste"))
force=ForceFromPotential("Harmonic")
params=LangevinParams(;n_steps = 10^5)

integrator_em=EM(force,1.0,0.001)
init_conds_em=initialize_initcond(integrator_em,init_conds_args)

integrator_bbk=BBK(force, 1.0, 1.0, 1.0, 1e-3)

# On peut tester aussi le calul des forces pour les différents moyen
const SUITE = BenchmarkGroup()

SUITE["init"] = BenchmarkGroup(["initialization"])

SUITE["init"]["overdamped"] = @benchmarkable generate_initial_conditions($integrator_em; params = $params,init_conds_args=$init_conds_args)
SUITE["init"]["underdamped"] = @benchmarkable generate_initial_conditions($integrator_bbk; params = $params,init_conds_args=$init_conds_args)

SUITE["integrator"] = BenchmarkGroup(["run"])

state_em=InitState(integrator_em, init_conds_em)
SUITE["integrator"]["run_EM"] = @benchmarkable run_trajectory!($state_em, $integrator_em; params = $params)

init_conds=initialize_initcond(integrator_bbk,init_conds_args)
state_bbk=InitState(integrator_bbk, init_conds)
SUITE["integrator"]["run_BBK"] = @benchmarkable run_trajectory!($state_bbk, $integrator_bbk; params = $params)

integrator_gjf=GJF(force, 1.0, 1.0, 1.0, 1e-3)
state_gjf=InitState(integrator_gjf, init_conds)
SUITE["integrator"]["run_GJF"] = @benchmarkable run_trajectory!($state_gjf, $integrator_gjf; params = $params)

integrator_aboba=ABOBA(force, 1.0, 1.0, 1.0, 1e-3)
state_aboba=InitState(integrator_aboba, init_conds)
SUITE["integrator"]["run_ABOBA"] = @benchmarkable run_trajectory!($state_aboba, $integrator_aboba; params = $params)

integrator_baoab=BAOAB(force, 1.0, 1.0, 1.0, 1e-3)
state_baoab=InitState(integrator_baoab, init_conds)
SUITE["integrator"]["run_BAOAB"] = @benchmarkable run_trajectory!($state_baoab, $integrator_baoab; params = $params)

integrator_verlet=Verlet(force, 1.0, 1e-3)
state_verlet=InitState(integrator_verlet, init_conds)
SUITE["integrator"]["run_Verlet"] = @benchmarkable run_trajectory!($state_verlet, $integrator_verlet; params = $params)


SUITE["forces"] = BenchmarkGroup(["run"])

force=ForceFromPotential("Harmonic")
integrator_em=EM(force,1.0,0.001)
init_conds_em=initialize_initcond(integrator_em,init_conds_args)
state_em=InitState(integrator_em, init_conds_em)

SUITE["forces"]["potential"] = @benchmarkable run_trajectory!($state_em, $integrator_em; params = $params)

coeffs=zeros(Float64,(1,2))
coeffs[1,:]=[0,-2]
force=ForceFromBasis("Taylor",coeffs)
integrator_em=EM(force,1.0,0.001)
init_conds_em=initialize_initcond(integrator_em,init_conds_args)
state_em=InitState(integrator_em, init_conds_em)

SUITE["forces"]["fromBasisTaylor"] = @benchmarkable run_trajectory!($state_em, $integrator_em; params = $params)
