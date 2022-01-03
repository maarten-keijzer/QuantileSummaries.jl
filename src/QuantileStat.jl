
# A Fast Algorithm for Approximate Quantiles in High Speed Data Streams
# Qi Zhang and Wei Wang

struct DataItem{T}
    value::T
    rmax::Int64
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
            rmax = xr.rmax + ys.rmax
        else
            rmax = xr.rmax + yt.rmax - 1
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

    count = list[end].rmax
    stepsize = count / sz

    newlist = QuantileSummary{T}(undef, sz+1)


    newlist[1] = list[1]
    i = 1
    next = 2;

    for rank = stepsize:stepsize:count
        while list[i].rmax < rank
            i += 1
        end

        newitem = DataItem{T}(
            list[i].value,
            list[i].rmax,
        )

        newlist[next] = newitem;
        next += 1
        i += 1
    end
    newlist
end

function percentile(qs::QuantileSummary, value)
    idx = searchsortedfirst(qs, DataItem(value, 0), by = x -> x.value)
    if idx > length(qs)
        idx = length(qs)
    end

    qs[idx].rmax / Float64(qs[end].rmax)
end

function qindex(qs::QuantileSummary, percentile::Number)
    idx = Int(round(percentile * length(qs))) + 1 # rank is one based
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


mutable struct QuantileStat{T} <: OnlineStat{T}
    eps::Float64

    current::Int64
    summary::QuantileSummary{T}
    workingset::FixedSizeSummaryBuilder{T}

    total::Int64
    countSoFar::Int64
    n::Int64

    cachedSummary::Union{Nothing, QuantileSummary{T}}

    QuantileStat{T}(eps::Float64 = 0.01) where {T} = begin
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
    QuantileStat(eps::Float64 = 0.01) = QuantileStat{Number}(eps)
end

function _fit!(qb::QuantileStat{T}, value) where {T}
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

targetsize(qb::QuantileStat) = Int64(floor(2/qb.eps))

function _merge!(qb1::QuantileStat{T}, qb2::QuantileStat{T}) where {T}
    
    # merge current summaries
    qb1.summary = merge(qb1.summary, qb2.summary)
    
    # merge the current working set of qb2
    qb1.summary = merge(qb1.summary, compress(merge(qb2.workingset), targetsize(qb1)))
    
    # Make sure nobs is still correct
    qb1.n += qb2.n

    # we keep qb1 in the same state as it was wrt exponentially increasing workingsets

    qb1
end

function summarize(qb::QuantileStat)::QuantileSummary
    if isnothing(qb.cachedSummary)
        qb.cachedSummary = merge(
            qb.summary,
            compress(merge(qb.workingset), Int64(floor(2 / qb.eps))),
        )

        # do a final compression to create a smaller summary that is still an ϵ-approximate summary
        qb.cachedSummary = compress(qb.cachedSummary, Int64(floor(2/qb.eps)))
    end
    qb.cachedSummary
end

ndata(qs::QuantileSummary) = length(qs)
ndata(qb::QuantileStat) = ndata(qb.summary) + ndata(qb.workingset)
ndata(fs::FixedSizeSummaryBuilder) = length(fs.values) + sum(ndata(qs) for qs in fs.levels)

value(qb::QuantileStat) = summarize(qb)

function qvalue(qb::QuantileStat, percentile::Number)
    qvalue(summarize(qb), percentile)
end

function qindex(qb::QuantileStat, percentile::Number)
    qindex(summarize(qb), percentile)
end

function percentile(qs::QuantileStat, value)
    percentile(summarize(qs), value)
end

function median(qb::QuantileStat)
    qvalue(qb, 0.5)
end
