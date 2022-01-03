
using Test
using Random
using QuantileSummaries


@testset "String" begin
    # string
    qb = QuantileBuilder{String}()

    for i = 1:1000
        fit!(qb, randstring(10))
    end

    summary = value(qb)
    @test true #passed without erroring
end

@testset "Correctness" begin
    
    eps = 0.01
    qb = QuantileBuilder(eps)

    n = 10_000
    r = randn(n)
    fit!(qb, r)
    sort!(r)

    qs = summarize(qb)

    # Test if the right percentile is identified (within eps)
    for percentile = eps:eps/2:1
        
        valueforpercentile = r[ Int(floor(percentile * n)) ]
        predictedpercentile = QuantileSummaries.percentile(qs, valueforpercentile)
        @test abs(predictedpercentile-percentile) < eps

    end

    for i = 1:1000
        percentile = rand() 
        idx = 1 + Int(floor(percentile * n))
        valueforpercentile = r[idx]
        predictedpercentile = QuantileSummaries.percentile(qs, valueforpercentile)
        @test abs(predictedpercentile-percentile) < eps
    end

    # check if computed value is within eps of expected
    # Percentiles & values are the same in this setup
    eps = 0.01
    qb = QuantileBuilder(eps)

    r = collect(0:0.0001:1)
    shuffle!(r)
    fit!(qb, r)
    qs = summarize(qb)

    for p in 0:0.001:1
        value = qvalue(qs, p)
        @test( abs(p-value) < eps )
    end
end

@testset "merge!" begin
    eps = 0.01
    qb1 = QuantileBuilder(eps)
    qb2 = QuantileBuilder(eps)

    n = 10_000
    r = randn(n)
    for value in r
        if rand() < 0.1
            fit!(qb1, value)
        else
            fit!(qb2, value)
        end
    end

    merge!(qb1, qb2)
    sort!(r)

    qs = summarize(qb1)
    for percentile = eps:eps/2:1
        valueforpercentile = r[ Int(floor(percentile * n)) ]
        predictedpercentile = QuantileSummaries.percentile(qs, valueforpercentile)
        @test abs(predictedpercentile-percentile) < eps
    end

end
