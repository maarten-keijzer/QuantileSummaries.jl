using QuantileSummaries
using Plots

qb = QuantileBuilder()

mem = []

function dorun(mem)
    for i = 1:500
        fit!(qb, randn(1_000_000))
        push!(mem, ndata(qb))
    end
end

plot(mem)
