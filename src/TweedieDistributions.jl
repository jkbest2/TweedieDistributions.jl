module TweedieDistributions

using Distributions
using Random

import Distributions: mean, var, succprob, failprob, insupport
import Base: rand, minimum, maximum

export Tweedie,
    CompoundPoissonGamma,
    rand,
    mean,
    var,
    succprob,
    failprob

abstract type AbstractTweedie{T} <: Distribution{Univariate,Continuous} end

"""
    Tweedie{T}
       μ::T
       p::T
       ϕ::T

Tweedie distribution. Typically used with 1 < ξ < 2 to give compound
gamma-Poisson. Has positive mass at zero, and is continuous for positive
Reals.
"""
struct Tweedie{T} <: AbstractTweedie{T}
    μ::T
    p::T
    ϕ::T
end

"""
    CompoundPoissonGamma{T} <: AbstractTweedie
        μ::T
        p::T
        ϕ::T

A Tweedie distribution with ``1 < p < 2`` is a compound Poisson-Gamma
distribution, making it more straightforward to draw random samples from.
"""
struct CompoundPoissonGamma{T} <: AbstractTweedie{T}
    μ::T
    p::T
    ϕ::T

    function CompoundPoissonGamma(μ::T, p::T, ϕ::T) where T
        μ ≥ 0 || throw("Mean must be nonnegative")
        ϕ > 0 || throw("Variance parameter must be positive")
        1 < p < 2 ||
            throw("CPG requires 1 < p < 2, otherwise use Tweedie.")
        new{T}(μ, p, ϕ)
    end
end

function CompoundPoissonGamma(Tw::Tweedie)
    CompoundPoissonGamma(Tw.μ, Tw.p, Tw.ϕ)
end

function CompoundPoissonGamma(μ, p, ϕ)
    μ_prom, p_prom, ϕ_prom = promote(μ, p, ϕ)
    CompoundPoissonGamma(μ_prom, p_prom, ϕ_prom)
end

@inline function _cpg_poissonmean(CPG::CompoundPoissonGamma)
    CPG.μ^(2 - CPG.p) / ((2 - CPG.p) * CPG.ϕ)
end

@inline function _cpg_gammashape(CPG::CompoundPoissonGamma)
    gam_shape = (2 - CPG.p) / (CPG.p - 1)
end

@inline function _cpg_gammarate(CPG::CompoundPoissonGamma)
    CPG.ϕ * (CPG.p - 1) * CPG.μ ^ (CPG.p - 1)
end

sampler(cpg::CompoundPoissonGamma) = CompoundPoissonGammaSampler(cpg)
rand(rng::AbstractRNG, cpg::CompoundPoissonGamma) = rand(rng, sampler(cpg))

# function rand(rng::Random.AbstractRNG, CPG::CompoundPoissonGamma{T}) where T
#     pois_mean = _cpg_poissonmean(CPG)
#     N = rand(rng, Poisson(pois_mean))
#     N == 0 && return zero(T)
#     gam_shape = N * _cpg_gammashape(CPG)
#     gam_rate = _cpg_gammarate(CPG)
# end

mean(Tw::AbstractTweedie) = Tw.μ
var(Tw::AbstractTweedie) = Tw.ϕ * Tw.μ ^ Tw.p

# rand(CPG::CompoundPoissonGamma) = rand(Random.GLOBAL_RNG, CPG)

function succprob(CPG::CompoundPoissonGamma{T}) where T
    N = _cpg_poissonmean(CPG)
    ccdf(Poisson(N), zero(T))
end

function failprob(CPG::CompoundPoissonGamma{T}) where T
    N = _cpg_poissonmean(CPG)
    cdf(Poisson(N), zero(T))
end

minimum(CPG::CompoundPoissonGamma{T}) where T = zero(T)
maximum(CPG::CompoundPoissonGamma) = Inf
# Need to restrict `x` to ::Real or dispatch is ambiguous
insupport(CPG::CompoundPoissonGamma, x::Real) = 0 ≤ x < Inf

struct CompoundPoissonGammaSampler{T} <: Sampleable{Univariate, Continuous}
    λ::T
    α::T
    θ::T
end

function CompoundPoissonGammaSampler(CPG::CompoundPoissonGamma{T}) where T
    λ = _cpg_poissonmean(CPG)
    α = _cpg_gammashape(CPG)
    θ = _cpg_gammarate(CPG)
    CompoundPoissonGammaSampler{T}(λ, α, θ)
end

function rand(rng::Random.AbstractRNG, cpgs::CompoundPoissonGammaSampler{T}) where T
    N = rand(rng, Poisson(cpgs.λ))
    N == 0 && return zero(T)
    rand(rng, Gamma(N * cpgs.α, cpgs.θ))
end
rand(cpgs::CompoundPoissonGammaSampler) = rand(Random.GLOBAL_RNG, cpgs)

end # module
