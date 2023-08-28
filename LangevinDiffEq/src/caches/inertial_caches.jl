
abstract type StochasticLangevinEqConstantCache <: StochasticDiffEqConstantCache end # Pourquoi faire ça, Si c'est pour avoir une seul function de check dans initialize!
abstract type StochasticLangevinEqMutableCache <: StochasticDiffEqMutableCache end


mutable struct GJConstantCache{uType,uMatrixType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::uMatrixType # constante de l'intégrateur
  sc₁::uMatrixType # constante de l'intégrateur
  d₁::uMatrixType # constante de l'intégrateur
  σ::uMatrixType # constante de l'intégrateur
end

@cache struct GJCache{uType,uMatrixType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  dutmp::uType
  k::uType
  gtmp::uType
  noise::rateNoiseType
  half::uEltypeNoUnits
  c₂::uMatrixType # constante de l'intégrateur
  sc₁::uMatrixType # constante de l'intégrateur
  d₁::uMatrixType # constante de l'intégrateur
  σ::uMatrixType # constante de l'intégrateur
  tmp::uTypeCombined
end

function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[2])
  c₂ = zero(noise_rate_prototype)
  sc₁ = zero(noise_rate_prototype)
  d₁ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)
  println("ConstantCache ",rate_prototype," ",noise_rate_prototype)
  GJConstantCache(k, uEltypeNoUnits(1//2), c₂, sc₁, d₁, σ)
end

function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  k = zero(rate_prototype.x[2])

  gtmp = zero(noise_rate_prototype)
  noise = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)
  c₂ = zero(noise_rate_prototype)
  sc₁ = zero(noise_rate_prototype)
  d₁ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)
  tmp = zero(u)
  println("Cache ",rate_prototype," ",noise_rate_prototype)
  GJCache(dutmp, k, gtmp, noise, half, c₂, sc₁, d₁, σ,tmp)
end
