using BenchmarkTools
using LangevinIntegrators
using Random

#Defining somes benchmark
# J'ai besoin de tester, l'évaluation d'une trajectoire pour les différents intégrateurs

# On peut tester aussi le calul des forces pour les différents moyen
const SUITE = BenchmarkGroup()

SUITE["trajectories"] = BenchmarkGroup(["string", "unicode"])
force=ForceFromPotential("Harmonic")
integrator=EM(force,1.0,0.001)
params=LangevinParams()
SUITE["trajectories"]["init_EM"] = @benchmarkable init_trajectory($integrator; params = $params)
state=init_trajectory(integrator; params = params)
SUITE["trajectories"]["run_EM"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=BBK(force, 1.0, 1.0, 1.0, 1e-3)
state=init_trajectory(integrator; params = params)
SUITE["trajectories"]["run_BBK"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=GJF(force, 1.0, 1.0, 1.0, 1e-3)
state=init_trajectory(integrator; params = params)
SUITE["trajectories"]["run_GJF"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=ABOBA(force, 1.0, 1.0, 1.0, 1e-3)
state=init_trajectory(integrator; params = params)
SUITE["trajectories"]["run_ABOBA"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=BAOAB(force, 1.0, 1.0, 1.0, 1e-3)
state=init_trajectory(integrator; params = params)
SUITE["trajectories"]["run_BAOAB"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

integrator=Verlet(force, 1.0, 1e-3)
state=init_trajectory(integrator; params = params)
SUITE["trajectories"]["run_Verlet"] = @benchmarkable run_trajectory!($state, $integrator; params = $params)

#
# SUITE["trigonometry"] = BenchmarkGroup(["math", "triangles"])
# SUITE["trigonometry"]["circular"] = BenchmarkGroup()
# for f in (sin, cos, tan)
#     for x in (0.0, pi)
#         SUITE["trigonometry"]["circular"][string(f), x] = @benchmarkable ($f)($x)
#     end
# end
#
# SUITE["trigonometry"]["hyperbolic"] = BenchmarkGroup()
# for f in (sin, cos, tan)
#     for x in (0.0, pi)
#         SUITE["trigonometry"]["hyperbolic"][string(f), x] = @benchmarkable ($f)($x)
#     end
# end
