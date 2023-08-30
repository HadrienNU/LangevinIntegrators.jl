
abstract type StochasticLangevinEqConstantCache <: StochasticDiffEqConstantCache end # Pourquoi faire ça, Si c'est pour avoir une seul function de check dans initialize!
abstract type StochasticLangevinEqMutableCache <: StochasticDiffEqMutableCache end


mutable struct GJConstantCache{uType, rateNoiseType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  type::String
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::rateNoiseType
  sc₁::rateNoiseType
  d₁::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits

end

@cache struct GJCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  type::String
  dutmp::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::rateNoiseType
  sc₁::rateNoiseType
  d₁::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
  tmp::uTypeCombined
end

mutable struct GJNonDiagonalNoiseConstantCache{uType,vMatrixType,rateNoiseType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  type::String
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::rateNoiseType
  sc₁::rateNoiseType
  d₁::vMatrixType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
end

@cache struct GJNonDiagonalNoiseCache{uType,vMatrixType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  type::String
  dutmp::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::rateNoiseType
  sc₁::rateNoiseType
  d₁::vMatrixType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
  tmp::uTypeCombined
end



function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  c₂ = zero(noise_rate_prototype)
  sc₁ = zero(noise_rate_prototype)
  d₁ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
      return GJConstantCache(alg.type,k, uEltypeNoUnits(1//2), c₂, sc₁, d₁, σ, sqrt(alg.kT))
  else
      return GJNonDiagonalNoiseConstantCache(alg.type,k, uEltypeNoUnits(1//2), c₂, sc₁, d₁, σ, sqrt(alg.kT))
  end

end

function alg_cache(alg::GJ,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  k = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)
  tmp = zero(u)

  noise = zero(rate_prototype.x[1])

  gtmp = zero(noise_rate_prototype)
  c₂ = zero(noise_rate_prototype)
  sc₁ = zero(noise_rate_prototype)
  d₁ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)

  if is_diagonal_noise(prob)
      return GJCache(alg.type,dutmp, k, gtmp, noise, half, c₂, sc₁, d₁, σ, sqrt(alg.kT),tmp)
  else
      return GJNonDiagonalNoiseCache(alg.type,dutmp, k, gtmp, noise, half, c₂, sc₁, d₁, σ, sqrt(alg.kT),tmp)
  end




  # @show ("Cache ",rate_prototype," ",noise_rate_prototype)

end
