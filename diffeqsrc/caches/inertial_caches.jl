
abstract type StochasticLangevinEqConstantCache <: StochasticDiffEqConstantCache end

mutable struct GJConstantCache{uType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::uEltypeNoUnits # constante de l'intégrateur
  sc₁::uEltypeNoUnits # constante de l'intégrateur
  d₁::uEltypeNoUnits # constante de l'intégrateur
  σ::uEltypeNoUnits # constante de l'intégrateur
end

abstract type StochasticLangevinEqMutableCache <: StochasticDiffEqMutableCache end

@cache struct GJCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticDiffEqMutableCache
  utmp::uType
  dutmp::uType
  k::uType
  gtmp::uType
  noise::rateNoiseType
  half::uEltypeNoUnits
  c₂::uEltypeNoUnits # constante de l'intégrateur
  sc₁::uEltypeNoUnits # constante de l'intégrateur
  d₁::uEltypeNoUnits # constante de l'intégrateur
  σ::uEltypeNoUnits # constante de l'intégrateur
end

function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  a= zero(rate_prototype.x[1])  # Il faudrait trouver plutot noise_rate_prototype ?
  c₂ = (1 - 0.5 * a) / (1 + 0.5 * a)
  sc₁ = sqrt((1 + c₂) / 2)
  d₁ = sqrt((1 - c₂) / a)
  σ = sqrt(2 * a)
  GJConstantCache(k, uEltypeNoUnits(1//2), uEltypeNoUnits(c₂), uEltypeNoUnits(sc₁), uEltypeNoUnits(d₁), uEltypeNoUnits(σ))
end

function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  gtmp = zero(rate_prototype.x[1])
  noise = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)
  c₂ = (1 - 0.5 * a) / (1 + 0.5 * a)
  sc₁ = sqrt((1 + c₂) / 2)
  d₁ = sqrt((1 - c₂) / a)
  σ = sqrt(2 * a)

  GJCache(utmp, dutmp, k, gtmp, noise, half, uEltypeNoUnits(c₂), uEltypeNoUnits(sc₁), uEltypeNoUnits(d₁), uEltypeNoUnits(σ))
end


# Fonction pour BAOAB tiré de SotchasticDiffEq.jl, à modifier

mutable struct BAOABConstantCache{uType,uEltypeNoUnits} <: StochasticDiffEqConstantCache
  k::uType  #
  half::uEltypeNoUnits  # =0.5
  c1::uEltypeNoUnits # constante de l'intégrateur
  c2::uEltypeNoUnits # constante de l'intégrateur
end
@cache struct BAOABCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticDiffEqMutableCache
  utmp::uType
  dutmp::uType
  k::uType
  gtmp::uType
  noise::rateNoiseType
  half::uEltypeNoUnits
  c1::uEltypeNoUnits
  c2::uEltypeNoUnits
  tmp::uTypeCombined
end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  c1 = exp(-alg.gamma*dt)
  c2 = sqrt(1 - alg.scale_noise*c1^2) # if scale_noise == false, c2 = 1
  BAOABConstantCache(k, uEltypeNoUnits(1//2), uEltypeNoUnits(c1), uEltypeNoUnits(c2))
end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  gtmp = zero(rate_prototype.x[1])
  noise = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)
  c1 = exp(-alg.gamma*dt)
  c2 = sqrt(1 - alg.scale_noise*c1^2) # if scale_noise == false, c2 = 1

  tmp = zero(u)

  BAOABCache(utmp, dutmp, k, gtmp, noise, half, uEltypeNoUnits(c1), uEltypeNoUnits(c2), tmp)
end
