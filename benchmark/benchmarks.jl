using BenchmarkTools
using LangevinIntegrators
# using Random

#Defining somes benchmark
# J'ai besoin de tester, l'évaluation d'une trajectoire pour les différents intégrateurs

benchmarkdir=abspath(dirname(@__FILE__ ))


init_conds_args=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"Cste"))
force=ForceFromPotential("Harmonic")
params=TrajsParams(;n_steps = 10^5)

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


params_full,init_conds_args_hidden=read_conf(joinpath(benchmarkdir, "hidden_comparison.ini"))

integrator_hidden=read_integrator_hidden_npz(joinpath(benchmarkdir, "coeffs_benchmark.npz");force=force)
init_conds_hidden=initialize_initcond(integrator_hidden,init_conds_args_hidden)
state_hidden=InitState(integrator_hidden, init_conds_hidden)

SUITE["integrator"]["run_EM_Hidden"] = @benchmarkable run_trajectory!($state_hidden, $integrator_hidden; params = $params)

integrator_aboba_hidden=read_integrator_hidden_npz(joinpath(benchmarkdir, "coeffs_benchmark.npz"); integrator_type ="ABOBA", force=force)
init_conds_hidden=initialize_initcond(integrator_aboba_hidden,init_conds_args_hidden)
state_aboba_hidden=InitState(integrator_aboba_hidden, init_conds_hidden)

SUITE["integrator"]["run_ABOBAHidden"] = @benchmarkable run_trajectory!($state_aboba_hidden, $integrator_aboba_hidden; params = $params)

integrator_hidden=read_integrator_hidden_npz(joinpath(benchmarkdir, "coeffs_benchmark.npz"))

SUITE["fullset"] = BenchmarkGroup(["run"])
SUITE["fullset"]["hidden"] = @benchmarkable run_trajectories($integrator_hidden; params = $params_full, init_conds_args=$init_conds_args_hidden)

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

force= ForceFromScipySplines(3,[0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549],[8.80178094 -0.25194201 1.65292369 -8.13424986 -10.49673538])
integrator_em=EM(force,1.0,0.001)
init_conds_em=initialize_initcond(integrator_em,init_conds_args)
state_em=InitState(integrator_em, init_conds_em)

SUITE["forces"]["fromScipySplines"] = @benchmarkable run_trajectory!($state_em, $integrator_em; params = $params)


force= ForceFromSplines(3,[0.163005, 0.163005, 0.163005, 0.163005, 0.622277, 1.081549, 1.081549, 1.081549, 1.081549],[8.80178094 -0.25194201 1.65292369 -8.13424986 -10.49673538])
integrator_em=EM(force,1.0,0.001)
init_conds_em=initialize_initcond(integrator_em,init_conds_args)
state_em=InitState(integrator_em, init_conds_em)

SUITE["forces"]["fromBSplinesKit"] = @benchmarkable run_trajectory!($state_em, $integrator_em; params = $params)
