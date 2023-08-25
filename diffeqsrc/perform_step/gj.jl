# Contenu actuel pour GJ

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

function initialize!(integrator, cache::GJConstantCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, du1, u1, p, t, integrator, cache)
  cache.k .= integrator.f.f1(du1,u1,p,t)

  γ=integrator.f.g(du1,u1,p,t)

  a = γ * dt
  if integrator.alg.type == "I"
      c₂ = (1 - 0.5 * a) / (1 + 0.5 * a)
  elseif integrator.alg.type == "II"
      c₂ = exp(-a)
  elseif integrator.alg.type == "III"
      c₂ = 1 - a
  elseif integrator.alg.type == "IV"
      c₂ = (sqrt(1 + 4 * a) - 1) / (2 * a)
  elseif integrator.alg.type == "V"
      c₂ = 1 / (1 + a)
  elseif integrator.alg.type == "VI"
      c₂ = 1 / (1 + 0.5 * a)^2
  else # Raise an error
      println("Unknown GJ type")
  end
  cache.sc₁ = sqrt((1 + c₂) / 2)
  cache.d₁ = sqrt((1 - c₂) / a)
  cache.σ = sqrt(2 * γ * Δt / β)
end

function initialize!(integrator, cache::GJCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, cache.k, du1, u1, p, t, integrator, cache)
  integrator.f.f1(cache.k,du1,u1,p,t)

  γ=integrator.f.g(du1,u1,p,t)

  a = γ * dt
  if integrator.alg.type == "I"
      c₂ = (1 - 0.5 * a) / (1 + 0.5 * a)
  elseif integrator.alg.type == "II"
      c₂ = exp(-a)
  elseif integrator.alg.type == "III"
      c₂ = 1 - a
  elseif integrator.alg.type == "IV"
      c₂ = (sqrt(1 + 4 * a) - 1) / (2 * a)
  elseif integrator.alg.type == "V"
      c₂ = 1 / (1 + a)
  elseif integrator.alg.type == "VI"
      c₂ = 1 / (1 + 0.5 * a)^2
  else # Raise an error
      println("Unknown GJ type")
  end
  cache.sc₁ = sqrt((1 + c₂) / 2)
  cache.d₁ = sqrt((1 - c₂) / a)
  cache.σ = sqrt(2 * γ * Δt / β)
end

@muladd function perform_step!(integrator,cache::GJConstantCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack k, half, c1, c2 = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # B
  du2 = du1 + half*dt*k

  # A
  u2 = u1 + half*dt*du2

  # O
  noise = integrator.g(u2,p,t+dt*half).*W.dW / sqdt
  du3 = c1*du2 + c2*noise

  # A
  u = u2 + half*dt*du3

  # B
  k .= f.f1(du3,u,p,t+dt)
  du = du3 + half*dt*k

  integrator.u = ArrayPartition((du, u))
end

@muladd function perform_step!(integrator,cache::GJCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack utmp, dutmp, k, gtmp, noise, half, c1, c2 = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # B
  @.. dutmp = du1 + half*dt*k

  # A
  @.. utmp = u1 + half*dt*dutmp

  # O
  integrator.g(gtmp,utmp,p,t+dt*half)
  @.. noise = gtmp*W.dW / sqdt
  @.. dutmp = c1*dutmp + c2*noise

  # A
  @.. u.x[2] = utmp + half*dt*dutmp

  # B
  f.f1(k,dutmp,u.x[2],p,t+dt)
  @.. u.x[1] = dutmp + half*dt*k
end
