function verify_f2(f, p, q, pa, t, integrator, ::StochasticLangevinEqConstantCache)
    res = f(p, q, pa, t)
    res != p && throwex(integrator)
end

function verify_f2(f, res, p, q, pa, t, integrator, ::StochasticLangevinEqMutableCache)
    f(res, p, q, pa, t)
    res != p && throwex(integrator)
end
function throwex(integrator)
  algn = typeof(integrator.alg)
  throw(ArgumentError("Algorithm $algn is not applicable if f2(p, q, t) != p"))
end

function initialize!(integrator, cache::StochasticLangevinEqConstantCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, du1, u1, p, t, integrator, cache)
  cache.k = integrator.f.f1(du1,u1,p,t)

  γ=integrator.g(u1,p,t)

  langevin_coefficients(γ * dt, cache)

end

function initialize!(integrator, cache::StochasticLangevinEqMutableCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, cache.k, du1, u1, p, t, integrator, cache)
  integrator.f.f1(cache.k,du1,u1,p,t)

  integrator.g(du1,u1,p,t)

  langevin_coefficients(du1 .* dt, cache)
end
