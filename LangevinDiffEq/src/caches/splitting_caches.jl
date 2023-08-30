mutable struct BAOABConstantCache{uType,rateNoiseType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
end

@cache struct BAOABCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  dutmp::uType
  utmp::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
  tmp::uTypeCombined
end

mutable struct BAOABNonDiagonalNoiseConstantCache{uType,rateNoiseType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
end

@cache struct BAOABNonDiagonalNoiseCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  dutmp::uType
  utmp::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
  tmp::uTypeCombined
end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  c₂ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
      return BAOABConstantCache(k, uEltypeNoUnits(1//2), c₂, σ, sqrt(alg.kT))
  else
      return BAOABNonDiagonalNoiseConstantCache(k, uEltypeNoUnits(1//2), c₂, σ, sqrt(alg.kT))
  end

end

function alg_cache(alg::BAOAB,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)
  tmp = zero(u)

  noise = zero(rate_prototype.x[1])

  gtmp = zero(noise_rate_prototype)
  c₂ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)

  if is_diagonal_noise(prob)
      return BAOABCache(dutmp, utmp , k, gtmp, noise, half, c₂, σ, sqrt(alg.kT),tmp)
  else
      return BAOABNonDiagonalNoiseCache(dutmp, utmp, k, gtmp, noise, half, c₂,σ, sqrt(alg.kT),tmp)
  end
  # println("Cache ",rate_prototype," ",noise_rate_prototype)
end



mutable struct ABOBAConstantCache{uType,rateNoiseType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
end

@cache struct ABOBACache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  dutmp::uType
  utmp::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
  tmp::uTypeCombined
end

mutable struct ABOBANonDiagonalNoiseConstantCache{uType,rateNoiseType,uEltypeNoUnits} <: StochasticLangevinEqConstantCache
  k::uType  #Force at previous timestep
  half::uEltypeNoUnits  # =0.5
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
end

@cache struct ABOBANonDiagonalNoiseCache{uType,uEltypeNoUnits,rateNoiseType,uTypeCombined} <: StochasticLangevinEqMutableCache
  dutmp::uType
  utmp::uType
  k::uType
  gtmp::rateNoiseType
  noise::uType
  half::uEltypeNoUnits
  c₂::rateNoiseType
  σ::rateNoiseType
  sqkT::uEltypeNoUnits
  tmp::uTypeCombined
end

function alg_cache(alg::ABOBA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{false}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  k = zero(rate_prototype.x[1])
  c₂ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)
  if is_diagonal_noise(prob)
      return ABOBAConstantCache(k, uEltypeNoUnits(1//2), c₂, σ, sqrt(alg.kT))
  else
      return ABOBANonDiagonalNoiseConstantCache(k, uEltypeNoUnits(1//2), c₂, σ, sqrt(alg.kT))
  end

end

function alg_cache(alg::ABOBA,prob,u,ΔW,ΔZ,p,rate_prototype,noise_rate_prototype,jump_rate_prototype,::Type{uEltypeNoUnits},::Type{uBottomEltypeNoUnits},::Type{tTypeNoUnits},uprev,f,t,dt,::Type{Val{true}}) where {uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits}
  dutmp = zero(u.x[1])
  utmp = zero(u.x[2])
  k = zero(rate_prototype.x[1])

  half = uEltypeNoUnits(1//2)
  tmp = zero(u)

  noise = zero(rate_prototype.x[1])

  gtmp = zero(noise_rate_prototype)
  c₂ = zero(noise_rate_prototype)
  σ = zero(noise_rate_prototype)

  if is_diagonal_noise(prob)
      return ABOBACache(dutmp, utmp , k, gtmp, noise, half, c₂, σ, sqrt(alg.kT),tmp)
  else
      return ABOBANonDiagonalNoiseCache(dutmp, utmp, k, gtmp, noise, half, c₂,σ, sqrt(alg.kT),tmp)
  end
  # println("Cache ",rate_prototype," ",noise_rate_prototype)
end
