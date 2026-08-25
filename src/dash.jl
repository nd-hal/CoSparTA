_dash_gamma_len(g) = g === nothing ? 0 : length(g)

function _dash_resolve_covariate_names(covariate_names, gamma_lens, supervised, K,
                                        intercept_flag, intercept_specified)
    if covariate_names === nothing
        return (fill(nothing, K), intercept_flag)
    end

    is_vector_form = covariate_names isa AbstractVector{<:AbstractString}
    is_list_form = !is_vector_form && covariate_names isa AbstractVector &&
        all(x -> x === nothing || x isa AbstractVector{<:AbstractString}, covariate_names)

    if !is_vector_form && !is_list_form
        error("covariate_names must be a vector of strings or a length-K list of per-component names.")
    end

    names = Vector{Union{Nothing,Vector{String}}}(fill(nothing, K))

    if is_list_form
        length(covariate_names) == K ||
            error("covariate_names supplied as a list must have length K=$K; got $(length(covariate_names)).")

        if !intercept_specified
            named_k = [k for k in 1:K if supervised[k] && covariate_names[k] !== nothing]
            if !isempty(named_k)
                matches_no_intercept = all(k -> length(covariate_names[k]) == gamma_lens[k] - 1, named_k)
                matches_with_intercept = all(k -> length(covariate_names[k]) == gamma_lens[k], named_k)
                if !matches_no_intercept && matches_with_intercept
                    intercept_flag = false
                    @warn "covariate_names lengths match gamma lengths exactly, implying no intercept; proceeding with intercept=false. Pass intercept=false explicitly to silence this warning."
                end
            end
        end

        for k in 1:K
            if supervised[k] && covariate_names[k] !== nothing
                expected = intercept_flag ? gamma_lens[k] - 1 : gamma_lens[k]
                length(covariate_names[k]) == expected ||
                    error("covariate_names for component $k must have length $expected; got $(length(covariate_names[k])).")
                names[k] = collect(String, covariate_names[k])
            end
        end

        return (names, intercept_flag)
    end

    supervised_lens = unique(gamma_lens[supervised])
    if length(supervised_lens) > 1
        error("covariate_names supplied as a single vector but components have different covariate counts (list them); supply a length-K list instead.")
    end
    isempty(supervised_lens) && return (names, intercept_flag)

    L = supervised_lens[1]

    if !intercept_specified && length(covariate_names) == L && length(covariate_names) != L - 1
        intercept_flag = false
        @warn "covariate_names lengths match gamma lengths exactly, implying no intercept; proceeding with intercept=false. Pass intercept=false explicitly to silence this warning."
    end

    expected = intercept_flag ? L - 1 : L
    length(covariate_names) == expected ||
        error("covariate_names must have length $expected; got $(length(covariate_names)).")

    cn = collect(String, covariate_names)
    for k in 1:K
        supervised[k] && (names[k] = cn)
    end

    return (names, intercept_flag)
end

function _dash_component_covariate_names(k, gamma_lens, intercept_flag, resolved_names)
    q = intercept_flag ? gamma_lens[k] - 1 : gamma_lens[k]
    q <= 0 && return String[]
    resolved_names[k] !== nothing && return resolved_names[k]
    return ["V$i" for i in 1:q]
end

function _dash_all_pairs_disjoint(name_sets)
    n = length(name_sets)
    for i in 1:n, j in (i+1):n
        isempty(intersect(name_sets[i], name_sets[j])) || return false
    end
    return true
end

function _dash_build_gamma_table(gamma_list, gamma_lens, supervised, resolved_names, K, intercept_flag)
    comp_names = [_dash_component_covariate_names(k, gamma_lens, intercept_flag, resolved_names) for k in 1:K]
    supervised_idx = findall(supervised)
    mode = _dash_all_pairs_disjoint([comp_names[k] for k in supervised_idx]) ? :pairs : :union

    rank = String[]
    note = String[]
    intercept_col = intercept_flag ? Union{Missing,Float64}[] : nothing

    if mode === :union
        gamma_col_names = String[]
        for k in supervised_idx, nm in comp_names[k]
            nm in gamma_col_names || push!(gamma_col_names, nm)
        end
        value_cols = Dict(nm => Union{Missing,Float64}[] for nm in gamma_col_names)

        for k in 1:K
            push!(rank, "R$k")
            if supervised[k]
                push!(note, "")
                covariate_vals = intercept_flag ? gamma_list[k][2:end] : gamma_list[k]
                intercept_flag && push!(intercept_col, gamma_list[k][1])
                row_vals = Dict{String,Any}(nm => missing for nm in gamma_col_names)
                for (i, nm) in enumerate(comp_names[k])
                    row_vals[nm] = covariate_vals[i]
                end
                for nm in gamma_col_names
                    push!(value_cols[nm], row_vals[nm])
                end
            else
                push!(note, "no covariates for this component")
                intercept_flag && push!(intercept_col, missing)
                for nm in gamma_col_names
                    push!(value_cols[nm], missing)
                end
            end
        end

        df = DataFrame(Rank=rank, Note=note)
        intercept_flag && (df.Intercept = intercept_col)
        for nm in gamma_col_names
            df[!, Symbol(nm)] = value_cols[nm]
        end
        return (df, :union)
    else
        covariates_col = String[]
        for k in 1:K
            push!(rank, "R$k")
            if supervised[k]
                push!(note, "")
                covariate_vals = intercept_flag ? gamma_list[k][2:end] : gamma_list[k]
                intercept_flag && push!(intercept_col, gamma_list[k][1])
                pieces = ["$(nm): $(round(covariate_vals[i]; digits=3))" for (i, nm) in enumerate(comp_names[k])]
                push!(covariates_col, join(pieces, ", "))
            else
                push!(note, "no covariates for this component")
                intercept_flag && push!(intercept_col, missing)
                push!(covariates_col, "")
            end
        end

        df = DataFrame(Rank=rank, Note=note)
        intercept_flag && (df.Intercept = intercept_col)
        df.Covariates = covariates_col
        return (df, :pairs)
    end
end

"""Resolve dash() inputs (fit or raw factors, labels, covariate names) into a static gamma table."""
function dash_data(; fit=nothing, Ef=nothing, Ew=nothing, lambda=nothing,
                    gamma_list=nothing, channel_names=nothing, time_labels=nothing,
                    channel_groups=nothing, covariate_names=nothing,
                    intercept::Union{Bool,Nothing}=nothing)
    intercept_specified = intercept !== nothing
    intercept_flag = intercept === nothing ? true : intercept

    raw_supplied = Ef !== nothing || Ew !== nothing || lambda !== nothing

    local K
    if raw_supplied
        (Ef !== nothing && Ew !== nothing && lambda !== nothing) ||
            error("When bypassing fit, all of Ef, Ew, and lambda must be supplied together.")
        Ef isa AbstractMatrix || error("Ef must be a matrix.")
        Ew isa AbstractMatrix || error("Ew must be a matrix.")
        K = size(Ef, 2)
        size(Ew, 2) == K ||
            error("Ew must have the same number of columns (K) as Ef; got $(size(Ew, 2)) vs $K.")
        length(lambda) == K ||
            error("lambda must have length K=$K; got $(length(lambda)).")
        if gamma_list !== nothing
            length(gamma_list) == K ||
                error("gamma_list must have length K=$K; got $(length(gamma_list)).")
        end
    elseif fit !== nothing
        nf = normalize_factors(fit)
        Ef = nf.Ef
        Ew = nf.Ew
        lambda = nf.lambda
        gamma_list = nf.gamma_list
        K = size(Ef, 2)
    else
        error("Either fit or all of Ef, Ew, and lambda must be provided.")
    end

    if channel_names !== nothing
        length(channel_names) == size(Ew, 1) ||
            error("channel_names must have length matching size(Ew, 1)=$(size(Ew, 1)); got $(length(channel_names)).")
    end
    if time_labels !== nothing
        length(time_labels) == size(Ef, 1) ||
            error("time_labels must have length matching size(Ef, 1)=$(size(Ef, 1)); got $(length(time_labels)).")
    end
    if channel_groups !== nothing
        length(channel_groups) == size(Ew, 1) ||
            error("channel_groups must have length matching size(Ew, 1)=$(size(Ew, 1)); got $(length(channel_groups)).")
    end

    gamma_lens = gamma_list === nothing ? zeros(Int, K) : [_dash_gamma_len(gamma_list[k]) for k in 1:K]
    supervised = gamma_lens .> 0

    resolved_names, intercept_flag = _dash_resolve_covariate_names(
        covariate_names, gamma_lens, supervised, K, intercept_flag, intercept_specified)

    any_supervised = any(supervised)

    gamma_table, gamma_table_mode = if any_supervised
        _dash_build_gamma_table(gamma_list, gamma_lens, supervised, resolved_names, K, intercept_flag)
    else
        (nothing, :none)
    end

    return (Ef=Ef, Ew=Ew, lambda=lambda, gamma_list=gamma_list,
            channel_names=channel_names, time_labels=time_labels,
            channel_groups=channel_groups, covariate_names=resolved_names,
            intercept=intercept_flag, K=K, any_supervised=any_supervised,
            gamma_table=gamma_table, gamma_table_mode=gamma_table_mode)
end
