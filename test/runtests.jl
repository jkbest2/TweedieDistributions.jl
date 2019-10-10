using TweedieDistributions
using Test

Tw1 = Tweedie(1.0, 1.5, 1.0)
Tw2 = Tweedie(0.2, 0.2, 0.2)
CPG1 = CompoundPoissonGamma(Tw1)
CPG2 = CompoundPoissonGamma(2.0, 1.9, 1.9)

cpg1_rand = rand(CPG1)
cpg2_randvec = zeros(1_000)
for i in eachindex(cpg2_randvec)
    cpg2_randvec[i] = rand(CPG2)
end

@testset "TweedieDistributions.jl" begin
    @test cpg1_rand ≥ 0
    @test any(cpg2_randvec .== 0)
end
