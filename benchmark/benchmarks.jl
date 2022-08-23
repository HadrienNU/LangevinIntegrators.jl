using BenchmarkTools
using LangevinIntegrators
# using Random

#Defining somes benchmark
# J'ai besoin de tester, l'évaluation d'une trajectoire pour les différents intégrateurs

init_conds_args=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"Cste"))
force=ForceFromPotential("Harmonic")
params=LangevinParams(;n_steps = 10^4)

# On peut tester aussi le calul des forces pour les différents moyen
const SUITE = BenchmarkGroup()

SUITE["trajectories"] = BenchmarkGroup(["run"])

integrator=EM(force,1.0,0.001)
init_conds=initialize_initcond(integrator,init_conds_args)

SUITE["trajectories"]["init_EM"] = @benchmarkable InitState($integrator,$params)
state=InitState(integrator, init_conds)
SUITE["trajectories"]["run_EM"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=BBK(force, 1.0, 1.0, 1.0, 1e-3)
init_conds=initialize_initcond(integrator,init_conds_args)
state=InitState(integrator, init_conds)
SUITE["trajectories"]["run_BBK"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=GJF(force, 1.0, 1.0, 1.0, 1e-3)
state=InitState(integrator, init_conds)
SUITE["trajectories"]["run_GJF"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=ABOBA(force, 1.0, 1.0, 1.0, 1e-3)
state=InitState(integrator, init_conds)
SUITE["trajectories"]["run_ABOBA"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=BAOAB(force, 1.0, 1.0, 1.0, 1e-3)
state=InitState(integrator, init_conds)
SUITE["trajectories"]["run_BAOAB"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=Verlet(force, 1.0, 1e-3)
state=InitState(integrator, init_conds)
SUITE["trajectories"]["run_Verlet"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)
