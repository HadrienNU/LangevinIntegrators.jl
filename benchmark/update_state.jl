using LangevinIntegrators

#Defining somes benchmark

init_conds_args=Dict("position"=> Dict("type"=>"Cste"),"velocity"=> Dict("type"=>"Cste"))
force=ForceFromPotential("Harmonic")
params=TrajsParams(;n_steps = 10^5)
β = 1.0
Δt = 1e-3

integrator_em=EM(force,β,0.001)
init_conds_em=initialize_initcond(integrator_em,init_conds_args)

integrator_bbk=BBK(force, β, 1.0, 1.0, Δt)


state_em=InitState(integrator_em, init_conds_em)
UpdateState!(state_em, integrator_em)

init_conds=initialize_initcond(integrator_bbk,init_conds_args)
state_bbk=InitState(integrator_bbk, init_conds)
UpdateState!(state_bbk, integrator_bbk)

integrator_vec=VEC(force, β, 1.0, 1.0, Δt)
state_vec=InitState(integrator_vec, init_conds)
UpdateState!(state_vec, integrator_vec)


integrator_gjf=GJF(force, β, 1.0, 1.0, Δt)
state_gjf=InitState(integrator_gjf, init_conds)
UpdateState!(state_gjf, integrator_gjf)

integrator_aboba=ABOBA(force, β, 1.0, 1.0, Δt)
state_aboba=InitState(integrator_aboba, init_conds)
UpdateState!(state_aboba, integrator_aboba)

integrator_baoab=BAOAB(force, β, 1.0, 1.0, Δt)
state_baoab=InitState(integrator_baoab, init_conds)
UpdateState!(state_baoab, integrator_baoab)

integrator_obabo=OBABO(force, β, 1.0, 1.0, Δt)
state_obabo=InitState(integrator_obabo, init_conds)
UpdateState!(state_obabo, integrator_obabo)


integrator_verlet=Verlet(force, 1.0, Δt)
state_verlet=InitState(integrator_verlet, init_conds)
UpdateState!(state_verlet, integrator_verlet)


params_full,init_conds_args_hidden=read_conf("hidden_comparison.ini")

integrator_hidden=read_integrator_hidden_npz( "coeffs_benchmark.npz";force=force)
init_conds_hidden=initialize_initcond(integrator_hidden,init_conds_args_hidden)
state_hidden=InitState(integrator_hidden, init_conds_hidden)

UpdateState!(state_hidden, integrator_hidden)

integrator_aboba_hidden=read_integrator_hidden_npz("coeffs_benchmark.npz"; integrator_type ="ABOBA", force=force)
init_conds_hidden=initialize_initcond(integrator_aboba_hidden,init_conds_args_hidden)
state_aboba_hidden=InitState(integrator_aboba_hidden, init_conds_hidden)

UpdateState!(state_aboba_hidden, integrator_aboba_hidden)


#Benchmark for Kernel integrator
kernel= exp.(-20.0*LinRange(0,500*Δt, 500))

integrator_kernelbbk=BBK_Kernel(force, β , kernel, 1.0, Δt, 1)
init_conds=initialize_initcond(integrator_kernelbbk,init_conds_args)
state_kernelbbk=InitState(integrator_kernelbbk, init_conds)
UpdateState!(state_kernelbbk, integrator_kernelbbk)

integrator_kernelgjf=GJF_Kernel(force, β , kernel, 1.0, Δt, 1)
init_conds=initialize_initcond(integrator_kernelgjf,init_conds_args)
state_kernelgjf=InitState(integrator_kernelgjf, init_conds)
UpdateState!(state_kernelgjf, integrator_kernelgjf)
