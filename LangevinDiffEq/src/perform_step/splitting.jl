@memoize function langevin_coefficients(a,cache::BAOABConstantCache)
    cache.c₂ = exp.(-a)
    cache.σ = sqrt.(1 .-cache.c₂.^2)
end


@memoize function langevin_coefficients(a,cache::BAOABCache)
    @. cache.c₂ = exp.(-a)
    @. cache.σ = sqrt(1 .-cache.c₂.^2 ) #/ β
end


@muladd function perform_step!(integrator,cache::BAOABConstantCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack half, c₂, σ, sqkT = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]


  # @show sc₁
  # B: vt+1/2
  du_mid =  du1 + half*dt .* cache.k

  # A : xt+1/2
  u_mid = u1 + half*dt*du_mid
  # Et en vrai il faudrait u_t+0.5 pour le rentrer dans langevin_coefficients
  langevin_coefficients(integrator.g(u_mid,p,t+dt*half) * dt,cache)

  noise = sqkT*σ.*W.dW /sqdt

  # O
  du_mid_O = c₂ .* du_mid .+ noise

  # A: xt+1
  u = u_mid .+ half * dt .*du_mid_O

  # B: vt+1
  cache.k = f.f1(du_mid,u,p,t+dt)
  du = du_mid_O + half * dt .*cache.k

  integrator.u = ArrayPartition((du, u))

end

@muladd function perform_step!(integrator,cache::BAOABCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack dutmp, utmp, k, gtmp, noise, half, c₂, σ, sqkT = cache

  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # B: vt+1/2
  @.. dutmp = du1 + half*dt*k

  # A: xt+1/2
  @.. utmp = u1 + half*dt*dutmp

  # O
  integrator.g(gtmp,utmp,p,t+dt*half)
  langevin_coefficients(gtmp * dt,cache)
  @.. noise = sqkT*σ.*W.dW/sqdt
  @.. dutmp = c₂ .* dutmp .+ noise

  # A: xt+1
  @.. u.x[2] = utmp + half*dt*dutmp

  # B: vt+1
  f.f1(k,dutmp,u.x[2],p,t+dt)
  @.. u.x[1] = dutmp + half*dt*k

end



@memoize function langevin_coefficients(a,cache::ABOBAConstantCache)
    cache.c₂ = exp.(-a)
    cache.σ = sqrt.(1 .-cache.c₂.^2)
end


@memoize function langevin_coefficients(a,cache::ABOBACache)
    @. cache.c₂ = exp.(-a)
    @. cache.σ = sqrt(1 .-cache.c₂.^2 ) #/ β
end

@muladd function perform_step!(integrator,cache::ABOBAConstantCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack half, c₂, σ, sqkT = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]


  # @show sc₁
  # A : xt+1/2
  u_mid = u1 + half*dt*du_mid

  cache.k = f.f1(du1,u_mid,p,t+half*dt)

  langevin_coefficients(integrator.g(u_mid,p,t+dt*half) * dt,cache)
  noise = sqkT*σ.*W.dW /sqdt
  # BOB: vt+1
  du = c₂ .* (du1 + half*dt .* cache.k) .+ noise .+ half * dt .*cache.k

  # A: xt+1
  u = u_mid .+ half * dt .*du

  integrator.u = ArrayPartition((du, u))

end

@muladd function perform_step!(integrator,cache::ABOBACache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack dutmp, utmp, k, gtmp, noise, half, c₂, σ, sqkT = cache

  du1 = uprev.x[1]
  u1 = uprev.x[2]

  # A: xt+1/2
  @.. utmp = u1 + half*dt*du1

  # BOB: vt+1/2
  f.f1(k,du1,utmp,p,t+dt)
  # B
  @.. dutmp = du1 + half*dt*k
  # O
  integrator.g(gtmp,utmp,p,t+dt*half)
  langevin_coefficients(gtmp * dt,cache)
  @.. noise = sqkT*σ.*W.dW/sqdt
  @.. dutmp = c₂ .* dutmp .+ noise
  # B: vt+1

  @.. u.x[1] = dutmp + half*dt*k

  # A: xt+1
  @.. u.x[2] = utmp + half*dt*dutmp

end
