
abstract type StochasticLangevinEqConstantCache <: StochasticDiffEqConstantCache end # Pourquoi faire ça, Si c'est pour avoir une seul function de check dans initialize!
abstract type StochasticLangevinEqMutableCache <: StochasticDiffEqMutableCache end


mutable struct GJConstantCache{uType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::uType # constante de l'intégrateur
  sc₁::uType # constante de l'intégrateur
  d₁::uType # constante de l'intégrateur
  σ::uType # constante de l'intégrateur
end

@cache struct GJCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  utmp::uType
  dutmp::uType
  k::uType
  gtmp::uType
  noise::rateNoiseType
  half::uEltypeNoUnits
  c₂::uType # constante de l'intégrateur
  sc₁::uType # constante de l'intégrateur
  d₁::uType # constante de l'intégrateur
  σ::uType # constante de l'intégrateur
  tmp::uTypeCombined
end

function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  a= zero(rate_prototype.x[1])  # Il faudrait trouver plutot noise_rate_prototype ?
  c₂ = (1 .- 0.5 * a) ./ (1 .+ 0.5 * a)
  sc₁ = sqrt.((1 .+ c₂) / 2)
  d₁ = sqrt.((1 .- c₂) ./ a)
  σ = sqrt.(2 .* a)

  GJConstantCache(k, uEltypeNoUnits(1//2), c₂, sc₁, d₁, σ)
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
  tmp = zero(u)
  GJCache(utmp, dutmp, k, gtmp, noise, half, c₂, sc₁, d₁, σ,tmp)
end
