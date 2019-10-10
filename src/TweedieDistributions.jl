module TweedieDistributions

using Distributions
using Random

import Base.rand

export Tweedie,
    CompoundPoissonGamma,
    rand

abstract type AbstractTweedie <: Distribution{Univariate,Continuous} end

"""
    Tweedie{T}
       μ::T
       p::T
       ϕ::T

Tweedie distribution. Typically used with 1 < ξ < 2 to give compound
gamma-Poisson. Has positive mass at zero, and is continuous for positive
Reals.
"""
struct Tweedie{T} <: AbstractTweedie
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
struct CompoundPoissonGamma{T} <: AbstractTweedie
    μ::T
    p::T
    ϕ::T

    function CompoundPoissonGamma(μ::T, p::T, ϕ::T) where T
        μ > 0 || throw("Mean must be positive")
        ϕ > 0 || throw("Variance parameter must be positive")
        1 < p < 2 ||
            throw("CPG requires 1 < p < 2, otherwise use more general Tweedie.")
        new{T}(μ, p, ϕ)
    end
end


function CompoundPoissonGamma(Tw::Tweedie)
    CompoundPoissonGamma(Tw.μ, Tw.p, Tw.ϕ)
end

function CompoundPoissonGamma(μ, p, ϕ)
    μ_prom, p_prom, ϕ_prom = promote(μ, p, ϕ)
end

function rand(rng::Random.AbstractRNG, CPG::CompoundPoissonGamma{T}) where T
    pois_mean = CPG.μ^(2 - CPG.p) / ((2 - CPG.p) * CPG.ϕ)
    N = rand(rng, Poisson(pois_mean))
    N == 0 && return zero(T)
    gam_shape = N * (2 - CPG.p) / (CPG.p - 1)
    gam_rate = CPG.ϕ * (CPG.p - 1) * CPG.μ ^ (CPG.p - 1)
    rand(rng, Gamma(gam_shape, gam_rate))
end

rand(CPG::CompoundPoissonGamma) = rand(Random.GLOBAL_RNG, CPG)

end # module
