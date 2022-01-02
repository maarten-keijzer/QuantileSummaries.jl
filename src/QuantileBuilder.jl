

struct DataItem{T}
    value::T
    rankupperbound::Int64
end

const QuantileSummary{T} = Array{DataItem{T},1}

function init(values::Array{T,1}) where {T}
    sort!(values)
    summaries = QuantileSummary{T}()
    for i = 1:length(values)
        item = DataItem{T}(values[i], i)
        push!(summaries, item)
    end
    summaries
end

function Base.merge(x::QuantileSummary{T}, y::QuantileSummary{T})::QuantileSummary{T} where {T}

    if length(x) == 0
        return y
    end
    if length(y) == 0
        return x
    end

    idx = 1;
    list = QuantileSummary{T}(undef, length(x) + length(y))

    xi = 1
    yi = 1

    xprevious = nothing
    yprevious = nothing

    while xi <= length(x) || yi <= length(y)
        xitem = xi <= length(x) ? x[xi] : nothing
        yitem = yi <= length(y) ? y[yi] : nothing

        # find the insertion point

        if isnothing(xitem)
            xr = yitem
            ys = xprevious
            yt = nothing

            yprevious = yitem
            yi += 1

        elseif isnothing(yitem)
            xr = xitem
            ys = yprevious
            yt = nothing

            xprevious = xitem
            xi += 1
        else
            if xitem.value < yitem.value
                xr = xitem
                ys = yprevious
                yt = yitem

                xprevious = xitem
                xi += 1
            else

                xr = yitem
                ys = xprevious
                yt = xitem

                yprevious = yitem
                yi += 1

            end
        end

        if isnothing(yt)
            rmax = xr.rankupperbound + ys.rankupperbound
        else
            rmax = xr.rankupperbound + yt.rankupperbound - 1
        end

        list[idx] = DataItem{T}(xr.value, rmax);
        idx += 1;
    end
    list
end

function compress(list::QuantileSummary{T}, sz::Int)::QuantileSummary{T} where {T}

    if length(list) <= sz
        return list
    end

    count = list[end].rankupperbound
    stepsize = count / sz

    newlist = QuantileSummary{T}(undef, sz+1)

    i = 1

    newlist[1] = list[i]
    idx = 2;
    #push!(newlist, list[i])

    for rank = stepsize:stepsize:count
        while list[i].rankupperbound < rank
            i += 1
        end

        newitem = DataItem{T}(
            list[i].value,
            list[i].rankupperbound,
        )

        #push!(newlist, newitem)
        newlist[idx] = newitem;
        idx += 1
        i += 1
    end
    newlist
end

function maxdistance(l1::QuantileSummary{T}, l2::QuantileSummary{T} ) where T

    rankupperbound1 = Float64(l1[end].rankupperbound)
    rankupperbound2 = Float64(l2[end].rankupperbound)

    dist = 0.0

    upper = 1

    for item in l1
        while upper <= length(l2) && l2[upper].value < item.value
            upper += 1
        end
        if upper > length(l2)
            upper = length(l2)
        end

        p0 = item.rankupperbound / rankupperbound1
        if item.value == l2[upper].value
            p1 = l2[upper].rankupperbound / rankupperbound2
            dist = max(dist, abs(p0-p1))
        else
            if upper == 1
                upper = 2
            end

            if true
                #interpolate
                α = (item.value - l2[upper-1].value) / (l2[upper].value - l2[upper-1].value)
                if item.value > l2[upper].value
                    α = 1
                end

                p1upper = l2[upper].rankupperbound / rankupperbound2
                p1lower = l2[upper-1].rankupperbound / rankupperbound2
                p1 = p1lower + α * (p1upper - p1lower)
                dist = max(dist, abs(p0-p1))
                #println("$dist $α $upper $p1lower $p1upper $(item.value)")
            else
                p11 = abs(p0 - l2[upper-1].rankupperbound / rankupperbound2)
                p12 = abs(p0 - l2[upper].rankupperbound / rankupperbound2)
                dist = max(dist, max(p11,p12))
                #dist = max(dist, )
            end
        end

    end

    dist
end

function ks_probability(l1::QuantileSummary{T}, l2::QuantileSummary{T} ) where T
    diff = max(maxdistance(l1,l2), maxdistance(l2,l1))

    n = l1[end].rankupperbound
    m = l2[end].rankupperbound

    min(1.0, 2exp(- diff * diff * 2m / (1+m/n)))

end

function percentile(qs::QuantileSummary, value)
    idx = searchsortedfirst(qs, DataItem(value, 0), by = x -> x.value)
    if idx > length(qs)
        idx = length(qs)
    end

    qs[idx].rankupperbound / Float64(qs[end].rankupperbound)
end

function qindex(qs::QuantileSummary, percentile::Number)
    rank = percentile * qs[end].rankupperbound + 1 # rank is one based
    idx = searchsortedfirst(
        qs,
        DataItem(0, Int64(floor(rank))),
        by = x -> x.rankupperbound,
    )
    if idx > length(qs)
        idx = length(qs)
    end
    idx
end

function qvalue(qs::QuantileSummary, percentile::Number)
    qs[qindex(qs, percentile)].value
end

mutable struct FixedSizeSummaryBuilder{T}
    values::Array{T,1}
    levels::Array{QuantileSummary{T},1}

    maxsize::Int64

    FixedSizeSummaryBuilder{T}(eps::Float64, N::Int) where {T} = begin
        maxsize = Int64(floor(log2(eps * N) / eps))
        if (maxsize < 1)
            error("eps $(eps) with N = $N leads to illegal size $(maxsize) ")
        end
        new{T}(Array{T, 1}(), Array{QuantileSummary{T}, 1}(), maxsize)
    end

end

function _fit!(fs::FixedSizeSummaryBuilder{T}, value) where {T}
    push!(fs.values, value)
    if length(fs.values) == fs.maxsize
        summary = init(fs.values)
        inserted = false
        for i = 1:length(fs.levels)
            other = fs.levels[i]
            if isempty(other)
                fs.levels[i] = summary
                inserted = true
                break
            else
                summary = merge(summary, other)
                summary = compress(summary, fs.maxsize)
                fs.levels[i] = QuantileSummary{T}()
            end
        end
        if (!inserted)
            push!(fs.levels, summary)
        end
        empty!(fs.values)

    end
end

function summarize(fs::FixedSizeSummaryBuilder)::QuantileSummary
    summary = init(fs.values)
    for qs in fs.levels
        summary = merge(summary, qs)
    end
    compress(summary, fs.maxsize)
end

function Base.merge(fs::FixedSizeSummaryBuilder)::QuantileSummary
    summary = init(fs.values)
    for qs in fs.levels
        summary = merge(summary, qs)
    end
    summary
end


mutable struct QuantileBuilder{T} <: OnlineStat{T}
    eps::Float64

    current::Int64
    summary::QuantileSummary{T}
    workingset::FixedSizeSummaryBuilder{T}

    total::Int64
    countSoFar::Int64
    n::Int64

    cachedSummary::Union{Nothing, QuantileSummary{T}}

    QuantileBuilder{T}(eps::Float64 = 0.01) where {T} = begin
        current = 3
        total = Int64(floor(((1 << current) / eps)))
        n = 0
        workingset = FixedSizeSummaryBuilder{T}(eps / 2, total)
        new{T}(
            eps,
            current,
            QuantileSummary{T}(),
            workingset,
            total,
            n,
            0,
            nothing,
        )
    end
    QuantileBuilder(eps::Float64 = 0.01) = QuantileBuilder{Number}(eps)
end

function _fit!(qb::QuantileBuilder{T}, value) where {T}
    qb.cachedSummary = nothing

    _fit!(qb.workingset, T(value))
    qb.n += 1
    qb.countSoFar += 1
    if qb.countSoFar == qb.total
        qb.current += 1
        qb.total = Int64(floor((1 << qb.current) / qb.eps))
        qb.countSoFar = 1
        qb.summary = merge(
            qb.summary,
            compress(merge(qb.workingset), Int64(floor(2 / qb.eps))),
        )
        qb.workingset = FixedSizeSummaryBuilder{T}(qb.eps / 2, qb.total)
    end
end

function summarize(qb::QuantileBuilder)::QuantileSummary
    if isnothing(qb.cachedSummary)
        qb.cachedSummary = merge(
            qb.summary,
            compress(merge(qb.workingset), Int64(floor(2 / qb.eps))),
        )
        qb.cachedSummary = compress(qb.cachedSummary, Int64(floor(2/qb.eps)))
    end
    qb.cachedSummary
end

ndata(qs::QuantileSummary) = length(qs)
ndata(qb::QuantileBuilder) = ndata(qb.summary) + ndata(qb.workingset)
ndata(fs::FixedSizeSummaryBuilder) = length(fs.values) + sum(ndata(qs) for qs in fs.levels)

value(qb::QuantileBuilder) = summarize(qb)

function qvalue(qs::QuantileBuilder, percentile::Number)
    qvalue(summarize(qs), percentile)
end

function qindex(qs::QuantileBuilder, percentile::Number)
    qindex(summarize(qs), percentile)
end

function percentile(qs::QuantileBuilder, value)
    percentile(summarize(qs), value)
end

function median(qb::QuantileBuilder)
    qvalue(qb, 0.5)
end
