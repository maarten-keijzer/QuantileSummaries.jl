module QuantileSummaries


export QuantileBuilder, QuantileSummary, DataItem,
    fit!,
    value,
    nobs,
    summarize,
    percentile,
    qindex,
    qvalue,
    ndata

import StatsBase: nobs, fit!, merge!
import OnlineStatsBase: value, OnlineStat, _fit!, _merge!


include("QuantileBuilder.jl")

end # module
