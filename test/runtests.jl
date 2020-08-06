using TweedieDistributions
using Distributions
using Test

import Random: seed!

seed!(12345)

Tw1 = Tweedie(1.0, 1.5, 1.0)
Tw2 = Tweedie(0.2, 0.2, 0.2)
CPG1 = CompoundPoissonGamma(Tw1)
CPG2 = CompoundPoissonGamma(1.0, 1.7, 1.9)

n = 1_000

cpg1_rand = rand(CPG1)
cpg2_randvec = zeros(n)
for i in eachindex(cpg2_randvec)
    cpg2_randvec[i] = rand(CPG2)
end
cpg2_ci = cdf(Normal(), 0.975) *
    sqrt(failprob(CPG2) * succprob(CPG2) / n)

## Note that simultaneously testing three 95% confidence intervals mean all
## tests are expected to pass ~85% of the time. If one of those three fail,
## re-run.
@testset "TweedieDistributions.jl" begin
    @test eltype(Tw1) == Float64
    @test eltype(CPG2) == Float64
    @test cpg1_rand ≥ 0
    @test any(cpg2_randvec .== 0)
    @test succprob(CPG1) == 1 - exp(-CPG1.μ^(2 - CPG1.p) / ((2 - CPG1.p) * CPG1.ϕ))
    @test abs(failprob(CPG2) - mean(cpg2_randvec .== 0)) < cpg2_ci
    @test abs(mean(CPG2) - mean(cpg2_randvec)) < sqrt(var(CPG2) / n)
    @test (n - 1) * var(cpg2_randvec) / quantile(Chisq(n - 1), 0.975) <
        var(CPG2) <
        (n - 1) * var(cpg2_randvec) / quantile(Chisq(n - 1), 0.025)
    @test minimum(CPG1) == 0
    @test maximum(CPG2) == Inf
    @test all(insupport.(CPG2, cpg2_randvec))
end
