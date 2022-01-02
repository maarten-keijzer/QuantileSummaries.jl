
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

for i = 0:0.01:1
    println(Int(floor(i*100)), "% ", qvalue(summary, i))
end


# string
qb = QuantileBuilder{String}()

for i = 1:1000
    fit!(qb, randstring(10))
end

summary = value(qb)



