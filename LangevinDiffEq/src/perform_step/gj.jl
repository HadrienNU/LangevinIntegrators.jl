@memoize function langevin_coefficients(a, cache::GJConstantCache)
    if cache.type == "I"
        cache.c₂ = (1 .- 0.5 * a) ./ (1 .+ 0.5 * a)
    elseif cache.type == "II"
        cache.c₂ = exp.(-a)
    elseif cache.type == "III"
        cache.c₂ = 1 .- a
    elseif cache.type == "IV"
        cache.c₂ = (sqrt.(1 .+ 4 * a) .- 1) ./ (2 * a)
    elseif cache.type == "V"
        cache.c₂ = 1 ./ (1 .+ a)
    elseif cache.type == "VI"
        cache.c₂ = 1 ./ (1 .+ 0.5 * a).^2
    else # Raise an error
        throw(ArgumentError("Unknown GJ type: $type"))
    end
    cache.sc₁ = sqrt.((1 .+ cache.c₂) ./ 2)
    cache.d₁ = sqrt.((1 .- cache.c₂) ./ a)
    cache.σ = sqrt.(2 .* a)
end


@memoize function langevin_coefficients(a, cache::GJCache)
    if cache.type == "I"
        @. cache.c₂ = (1 - 0.5 * a) / (1 + 0.5 * a)
    elseif cache.type == "II"
        @. cache.c₂ = exp(-a)
    elseif cache.type == "III"
        @. cache.c₂ = 1 - a
    elseif cache.type == "IV"
        @. cache.c₂ = (sqrt(1 + 4 * a) - 1) / (2 * a)
    elseif cache.type == "V"
        @. cache.c₂ = 1 / (1 + a)
    elseif cache.type == "VI"
        @. cache.c₂ = 1 / (1 + 0.5 * a)^2
    else # Raise an error
        throw(ArgumentError("Unknown GJ type: $type"))
    end
    @. cache.sc₁ = sqrt((1 + cache.c₂) / 2)
    @. cache.d₁ = sqrt((1 - cache.c₂) / a)
    @. cache.σ = sqrt(2 * a )
end

@muladd function perform_step!(integrator,cache::GJConstantCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack half, c₂, sc₁, d₁, σ, sqkT = cache
  du1 = uprev.x[1]
  u1 = uprev.x[2]

  langevin_coefficients(integrator.g(u1,p,t) * dt, cache)


  # @show sc₁
  noise = sqkT*σ.*W.dW/sqdt
  # vt+1/2
  du_mid =  sc₁ .* du1 + d₁ .* half*dt .* cache.k +  half .* d₁ .* noise

  # xt+1
  u = u1 .+ dt .*d₁ .*  du_mid

  # vt+1
  cache.k = f.f1(du_mid,u,p,t+dt)
  du = ( c₂ .* du_mid + half *  dt .* d₁ .* cache.k  + half .* d₁ .* noise ) ./ sc₁

  integrator.u = ArrayPartition((du, u))

end

@muladd function perform_step!(integrator,cache::GJCache,f=integrator.f)
  @unpack t,dt,sqdt,uprev,u,p,W = integrator
  @unpack dutmp, k, gtmp, noise, half, c₂, sc₁, d₁, σ, sqkT = cache

  du1 = uprev.x[1]
  u1 = uprev.x[2]

  integrator.g(gtmp,u1,p,t)

  langevin_coefficients(gtmp * dt, cache)

  # @show sc₁

  @.. noise = sqkT*σ*W.dW/sqdt

  # vt+1/2
  @.. dutmp =  sc₁ * du1 + d₁ * half*dt * k +  half * d₁ * noise

  # xt+1
  @.. u.x[2] = u1 + d₁ * dt * dutmp


  # vt+1
  k = f.f1(k,dutmp,u.x[2],p,t+dt)
  @.. u.x[1] = ( c₂ * dutmp + half *  dt * d₁ * k  + half * d₁ * noise ) / sc₁

end
