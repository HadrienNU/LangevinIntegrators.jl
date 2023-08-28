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
  cache.k = integrator.f.f1(du1,u1,p,t)

  γ=integrator.g(u1,p,t)

  a = γ * dt
  if integrator.alg.type == "I"
      c₂ = (1 .- 0.5 * a) ./ (1 .+ 0.5 * a)
  elseif integrator.alg.type == "II"
      c₂ = exp.(-a)
  elseif integrator.alg.type == "III"
      c₂ = 1 .- a
  elseif integrator.alg.type == "IV"
      c₂ = (sqrt.(1 .+ 4 * a) .- 1) ./ (2 * a)
  elseif integrator.alg.type == "V"
      c₂ = 1 ./ (1 .+ a)
  elseif integrator.alg.type == "VI"
      c₂ = 1 ./ (1 .+ 0.5 * a).^2
  else # Raise an error
      println("Unknown GJ type")
  end
  @. cache.sc₁ = sqrt((1 + c₂) / 2)
  @. cache.d₁ = sqrt((1 - c₂) / a)
  @. cache.σ = sqrt(2 * γ * dt)
end

function initialize!(integrator, cache::GJCache)
  @unpack t,dt,uprev,u,p,W = integrator
  du1 = integrator.uprev.x[1]
  u1 = integrator.uprev.x[2]

  verify_f2(integrator.f.f2, cache.k, du1, u1, p, t, integrator, cache)
  integrator.f.f1(cache.k,du1,u1,p,t)

  γ=integrator.g(du1,u1,p,t)

  a = γ .* dt
  if integrator.alg.type == "I"
      c₂ = (1 .- 0.5 * a) ./ (1 .+ 0.5 * a)
  elseif integrator.alg.type == "II"
      c₂ = exp.(-a)
  elseif integrator.alg.type == "III"
      c₂ = 1 .- a
  elseif integrator.alg.type == "IV"
      c₂ = (sqrt.(1 .+ 4 * a) .- 1) ./ (2 * a)
  elseif integrator.alg.type == "V"
      c₂ = 1 ./ (1 .+ a)
  elseif integrator.alg.type == "VI"
      c₂ = 1 ./ (1 .+ 0.5 * a).^2
  else # Raise an error
      println("Unknown GJ type")
  end
  @. cache.sc₁ = sqrt((1 + c₂) / 2)
  @. cache.d₁ = sqrt((1 - c₂) / a)
  @. cache.σ = sqrt(2 * γ * dt ) #/ β
end


@muladd function perform_step!(integrator,cache::GJConstantCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack half, c₂, sc₁, d₁, σ = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  noise = σ*W.dW / sqdt  # Vérifier le bruit notamment du point de vue du pas de temps
  # vt+1/2
  du_mid =  sc₁ * du1 + d₁ * half*dt * cache.k +  half * d₁ * noise

  # xt+1
  u = u1 + d₁ * dt * du_mid

  # vt+1
  cache.k = f.f1(du_mid,u,p,t+dt)
  du = ( c₂ * du_mid + half *  dt * d₁ * cache.k  + half * d₁ * noise ) / sc₁

  integrator.u = ArrayPartition((du, u))

end

@muladd function perform_step!(integrator,cache::GJCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack dutmp, k, gtmp, noise, half, c₂, sc₁, d₁, σ = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  integrator.g(gtmp,u1,p,t+dt*half)
  @.. noise = σ*W.dW / sqdt

  # vt+1/2
  @.. dutmp =  sc₁ * du1 + d₁ * half*dt * k +  half * d₁ * noise

  # xt+1
  @.. u.x[2] = u1 + d₁ * dt * dutmp


  # vt+1
  k = f.f1(k,dutmp,u.x[2],p,t+dt)
  @.. u.x[1] = ( c₂ * dutmp + half *  dt * d₁ * k  + half * d₁ * noise ) / sc₁

end
