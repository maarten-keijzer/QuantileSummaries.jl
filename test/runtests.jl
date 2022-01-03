
using Test
using Random
using QuantileSummaries


@test true

qb = QuantileBuilder()

fit!(qb, randn(100_000))
@test nobs(qb) == 100_000

summary = value(qb)

qvalue(summary, 0.5)
qindex(summary, 500)

for i = 0:0.05:1
    println(Int(floor(i*100)), "% ", qvalue(summary, i))
end


# string
qb = QuantileBuilder{String}()

for i = 1:1000
    fit!(qb, randstring(10))
end

summary = value(qb)

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
        
        truevalue = r[ Int(floor(percentile * n)) ]
        predictedpercentile = QuantileSummaries.percentile(qs, truevalue)
        @test abs(predictedpercentile-percentile) < eps

    end

    for i = 1:1000
        percentile = rand() 
        idx = 1 + Int(floor(percentile * n))
        truevalue = r[idx]
        predictedpercentile = QuantileSummaries.percentile(qs, truevalue)
        @test abs(predictedpercentile-percentile) < eps

    end

end

@testset "Merge" begin
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
        
        truevalue = r[ Int(floor(percentile * n)) ]
        predictedpercentile = QuantileSummaries.percentile(qs, truevalue)
        @test abs(predictedpercentile-percentile) < eps

    end

end
