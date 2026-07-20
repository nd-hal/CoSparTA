function _unique_preserve(values)
    out = eltype(values)[]
    seen = Set{eltype(values)}()
    for x in values
        if !(x in seen)
            push!(out, x); push!(seen, x)
        end
    end
    out
end

function _bin_values(values, time_bins)
    temporal = eltype(values) <: Union{Date,DateTime}
    numeric = temporal ? Float64.(Dates.value.(values)) : Float64.(values)
    label_origin = temporal ? minimum(values) : nothing
    if time_bins isa Integer
        time_bins > 0 || throw(ArgumentError("time_bins must be positive"))
        lo, hi = extrema(numeric)
        if lo == hi
            return ones(Int, length(values)), ["[$lo, $hi]"]
        end
        breaks = collect(range(lo, hi; length=time_bins + 1))
    else
        raw_breaks = collect(time_bins)
        breaks = temporal ? Float64.(Dates.value.(raw_breaks)) : Float64.(raw_breaks)
        temporal && (label_origin = first(raw_breaks))
        length(breaks) >= 2 || throw(ArgumentError("explicit time-bin breaks need at least two values"))
        issorted(breaks) || throw(ArgumentError("time-bin breaks must be sorted"))
        minimum(numeric) >= first(breaks) && maximum(numeric) <= last(breaks) ||
            throw(ArgumentError("time-bin breaks must cover all values"))
    end
    nbins = length(breaks) - 1
    indices = [clamp(searchsortedlast(breaks, x), 1, nbins) for x in numeric]
    if temporal
        T = eltype(values)
        offsets = round.(Int, breaks .- first(breaks))
        restored = T <: DateTime ? label_origin .+ Millisecond.(offsets) : label_origin .+ Day.(offsets)
        formatted = string.(restored)
        labels = ["[$(formatted[i]), $(formatted[i + 1]))" for i in 1:nbins]
        labels[end] = "[$(formatted[end - 1]), $(formatted[end])]"
    else
        labels = ["[$(breaks[i]), $(breaks[i + 1]))" for i in 1:nbins]
        labels[end] = "[$(breaks[end - 1]), $(breaks[end])]"
    end
    return indices, labels
end

"""
Build a dense three-way count tensor from any Tables.jl-compatible long table.

`row`, `col`, `slice`, and `value` are symbols or strings. If `value=nothing`,
rows are counted as events. `time_bins` may be an integer or explicit numeric
break points.
"""
function build_tensor(data; row, col, slice, value=nothing, time_bins=nothing,
                      fill=0, row_levels=nothing, col_levels=nothing,
                      slice_levels=nothing, verbose=true)
    Tables.istable(data) || throw(ArgumentError("data must implement Tables.jl"))
    columns = Tables.columns(data)
    getcol(name) = collect(Tables.getcolumn(columns, Symbol(name)))
    rv = getcol(row); cv = getcol(col); sv = getcol(slice)
    N = length(rv)
    length(cv) == N == length(sv) || throw(DimensionMismatch("table columns have different lengths"))
    vv = value === nothing ? ones(Int, N) : getcol(value)
    length(vv) == N || throw(DimensionMismatch("value column has the wrong length"))
    (any(ismissing, rv) || any(ismissing, cv) || any(ismissing, sv) || any(ismissing, vv)) &&
        throw(ArgumentError("tensor index and value columns cannot contain missing"))

    if time_bins !== nothing
        cv, generated_col_levels = _bin_values(cv, time_bins)
        col_levels === nothing && (col_levels = collect(1:length(generated_col_levels)))
        output_col_labels = generated_col_levels
    else
        output_col_labels = nothing
    end

    rlevels = row_levels === nothing ? _unique_preserve(rv) : collect(row_levels)
    slevels = slice_levels === nothing ? _unique_preserve(sv) : collect(slice_levels)
    clevels = if col_levels !== nothing
        collect(col_levels)
    else
        sort(unique(cv))
    end
    rmap = Dict(x => i for (i, x) in enumerate(rlevels))
    cmap = Dict(x => i for (i, x) in enumerate(clevels))
    smap = Dict(x => i for (i, x) in enumerate(slevels))
    aggregates = Dict{NTuple{3,Int},Float64}()
    for z in 1:N
        haskey(rmap, rv[z]) || throw(ArgumentError("row value $(rv[z]) is absent from row_levels"))
        haskey(cmap, cv[z]) || throw(ArgumentError("column value $(cv[z]) is absent from col_levels"))
        haskey(smap, sv[z]) || throw(ArgumentError("slice value $(sv[z]) is absent from slice_levels"))
        key = (rmap[rv[z]], cmap[cv[z]], smap[sv[z]])
        aggregates[key] = get(aggregates, key, 0.0) + Float64(vv[z])
    end
    all(isinteger, values(aggregates)) || throw(ArgumentError("aggregated tensor counts must be integers"))
    X = Base.fill(Int(fill), length(rlevels), length(clevels), length(slevels))
    for (idx, count_value) in aggregates
        X[idx...] = Int(round(count_value))
    end
    verbose && @info "Tensor built" dimensions=size(X) percent_zeros=100 * count(iszero, X) / length(X)
    labels2 = output_col_labels === nothing ? clevels : output_col_labels
    return TensorData(X, collect(rlevels), collect(labels2), collect(slevels))
end

build_tensor(data, row, col, slice; kwargs...) =
    build_tensor(data; row=row, col=col, slice=slice, kwargs...)
