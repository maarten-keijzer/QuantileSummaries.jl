module QuantileSummaries


export QuantileStat, QuantileSummary, DataItem,
    fit!,
    value,
    nobs,
    summarize,
    percentile,
    qindex,
    qvalue,
    ndata

import StatsBase: nobs, fit!, merge!, value
import OnlineStatsBase: value, OnlineStat, _fit!, _merge!

include("QuantileStat.jl")

end # module
