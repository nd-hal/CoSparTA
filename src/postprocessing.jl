"""Retrieve a factor matrix, component weights, or covariate coefficients."""
function get_loadings(result::CoSparTAFit, what="U1"; normalized=true)
    key = uppercase(String(what))
    if key == "U1"
        return normalized ? result.U1_normed.mean : result.U1.mean
    elseif key == "U2"
        return normalized ? result.U2_normed.mean : result.U2.mean
    elseif key == "U3"
        return normalized ? result.U3_normed.mean : result.U3.mean
    elseif key in ("WEIGHT", "WEIGHTS", "LAMBDA")
        return normalized ? result.weights : result.weights_raw
    elseif key == "GAMMA"
        priors = normalized ? result.priors1[result.order] : result.priors1
        return [p === nothing ? nothing : p.gamma for p in priors]
    end
    throw(ArgumentError("what should be one of U1, U2, U3, weight, or gamma"))
end

"""Return unit-norm factors ordered by descending component weight."""
normalize_factors(result::CoSparTAFit) =
    (El=result.U1_normed.mean, Ef=result.U2_normed.mean, Ew=result.U3_normed.mean,
     lambda=result.weights, order=result.order, lambda_order=result.order,
     gamma_list=get_loadings(result, "gamma"))

"""Project one or more `p × w` observations onto learned mode-2/3 factors."""
function project_tensor(Xnew::AbstractArray; fit=nothing, normalize=true,
                        Ef=nothing, Ew=nothing, lambda=nothing)
    single = ndims(Xnew) == 2
    ndims(Xnew) in (2, 3) || throw(DimensionMismatch("Xnew must be p×w or n×p×w"))
    X3 = single ? reshape(Xnew, 1, size(Xnew, 1), size(Xnew, 2)) : Xnew
    if Ef === nothing || Ew === nothing
        fit === nothing && throw(ArgumentError("provide fit or Ef, Ew, and lambda"))
        Ef = normalize ? fit.U2_normed.mean : fit.U2.mean
        Ew = normalize ? fit.U3_normed.mean : fit.U3.mean
        lambda = normalize ? fit.weights : nothing
    end
    size(X3, 2) == size(Ef, 1) || throw(DimensionMismatch("time dimension does not match Ef"))
    size(X3, 3) == size(Ew, 1) || throw(DimensionMismatch("channel dimension does not match Ew"))
    size(Ef, 2) == size(Ew, 2) || throw(DimensionMismatch("Ef and Ew need the same number of columns"))
    Xmat = reshape(Float64.(X3), size(X3, 1), :)
    projection = zeros(size(X3, 1), size(Ef, 2))
    for k in axes(Ef, 2)
        kernel = vec(Ef[:, k] * transpose(Ew[:, k]))
        projection[:, k] .= Xmat * kernel
    end
    lambda !== nothing && (projection .*= reshape(Float64.(lambda), 1, :))
    return single ? vec(projection) : projection
end

project_tensor(Xnew, fit::CoSparTAFit; kwargs...) = project_tensor(Xnew; fit=fit, kwargs...)

"""Reconstruct the fitted Poisson mean tensor."""
function reconstruct_tensor(result::CoSparTAFit; normalized=true)
    L = normalized ? result.U1_normed.mean : result.U1.mean
    F = normalized ? result.U2_normed.mean : result.U2.mean
    W = normalized ? result.U3_normed.mean : result.U3.mean
    weights = normalized ? result.weights : ones(size(L, 2))
    Xhat = zeros(size(L, 1), size(F, 1), size(W, 1))
    for k in axes(L, 2), m in axes(W, 1), j in axes(F, 1), i in axes(L, 1)
        Xhat[i, j, m] += weights[k] * L[i, k] * F[j, k] * W[m, k]
    end
    Xhat
end

"""
Native nonnegative Poisson-CP warm start.

This replaces the R port's Python/`pyCP_APR` bridge and returns `(L, F, W)`.
"""
function init_cpapr(X::AbstractArray, K::Integer; n_iters=150, random_state=42,
                    rng=MersenneTwister(random_state))
    Xf, mask = _clean_tensor(X, nothing)
    n, p, w = size(Xf)
    L = rand(rng, n, K) .+ 0.1; F = rand(rng, p, K) .+ 0.1; W = rand(rng, w, K) .+ 0.1
    I, J, M, values = _coordinates(Xf, mask)
    for _ in 1:n_iters
        alpha = calc_qz_sparse(I, J, M, log.(L .+ _TINY), log.(F .+ _TINY), log.(W .+ _TINY))
        for k in 1:K
            Ez = calc_EZ_3d(I, J, M, values, @view(alpha[:, k]), n, p, w)
            L[:, k] .= max.(Ez.rs ./ max(sum(F[:, k]) * sum(W[:, k]), _TINY), _TINY)
            F[:, k] .= max.(Ez.cs ./ max(sum(L[:, k]) * sum(W[:, k]), _TINY), _TINY)
            W[:, k] .= max.(Ez.zs ./ max(sum(L[:, k]) * sum(F[:, k]), _TINY), _TINY)
        end
    end
    return (L, F, W)
end

function _ols_summary(y, X)
    A = hcat(ones(length(y)), Matrix{Float64}(X))
    beta = pinv(A) * y
    residuals = y - A * beta
    df = max(length(y) - rank(A), 1)
    sigma2 = sum(abs2, residuals) / df
    covariance = sigma2 .* pinv(A' * A)
    se = sqrt.(max.(diag(covariance), 0))
    tstat = beta ./ max.(se, _TINY)
    pvalue = 2 .* ccdf.(TDist(df), abs.(tstat))
    return (coef=beta, se=se, t=tstat, pvalue=pvalue, residual_df=df)
end

"""Screen covariates by OLS on active log-loadings for each component."""
function select_covariates(; K::Integer, covariate_data, fit=nothing, El=nothing,
                           X=nothing, alpha=0.05, covariate_names=nothing,
                           verbose=true, fit_kwargs...)
    fit_unsup = fit
    if El === nothing
        if fit === nothing
            X === nothing && throw(ArgumentError("provide El, fit, or X"))
            fit_unsup = CoSparTA.fit(X, K; Xcov=nothing, fit_kwargs...)
        end
        El = fit_unsup.U1_normed.mean
    end
    size(El, 2) == K || throw(DimensionMismatch("El must have K columns"))
    Xcov = Matrix{Float64}(covariate_data)
    size(Xcov, 1) == size(El, 1) || throw(DimensionMismatch("covariates and El need the same rows"))
    names = covariate_names === nothing ? ["V$i" for i in axes(Xcov, 2)] : String.(covariate_names)
    selected = Vector{Vector{Int}}(undef, K); summaries = Vector{Any}(undef, K)
    for k in 1:K
        active = findall(>(0), @view El[:, k])
        length(active) > size(Xcov, 2) + 1 || throw(ArgumentError("too few active rows for component $k"))
        summary = _ols_summary(log.(@view(El[active, k])), Xcov[active, :])
        selected[k] = findall(<(alpha), summary.pvalue[2:end])
        summaries[k] = merge(summary, (selected_names=names[selected[k]],))
        verbose && @info "Covariate screen" component=k selected=names[selected[k]]
    end
    return (selected=selected, summaries=summaries, fit_unsupervised=fit_unsup)
end

function _assignment_max(S)
    K = size(S, 1)
    if K > 20
        @warn "K > 20: using greedy factor assignment"
        used = falses(K); perm = zeros(Int, K)
        for i in sortperm(vec(maximum(S; dims=2)); rev=true)
            scores = copy(@view S[i, :]); scores[used] .= -Inf
            perm[i] = argmax(scores); used[perm[i]] = true
        end
        return perm
    end
    states = 1 << K
    dp = fill(-Inf, states); parent_mask = zeros(Int, states); parent_col = zeros(Int, states)
    dp[1] = 0.0
    for mask in 0:(states - 1)
        row = count_ones(mask) + 1
        row > K && continue
        for col in 1:K
            bit = 1 << (col - 1)
            mask & bit != 0 && continue
            next = mask | bit
            score = dp[mask + 1] + S[row, col]
            if score > dp[next + 1]
                dp[next + 1] = score; parent_mask[next + 1] = mask; parent_col[next + 1] = col
            end
        end
    end
    perm = zeros(Int, K); mask = states - 1
    for row in K:-1:1
        perm[row] = parent_col[mask + 1]; mask = parent_mask[mask + 1]
    end
    perm
end

"""Optimally align estimated factor columns to reference factors."""
function match_factors(ref, est; absolute_value=true)
    length(ref) == length(est) || throw(DimensionMismatch("ref and est must have the same modes"))
    K = size(first(ref), 2); S = ones(K, K)
    for (A, B) in zip(ref, est)
        size(A) == size(B) || throw(DimensionMismatch("corresponding factor matrices must match"))
        size(A, 2) == K || throw(DimensionMismatch("all modes must have K columns"))
        An = Matrix(A) ./ reshape(max.([norm(A[:, k]) for k in 1:K], _TINY), 1, :)
        Bn = Matrix(B) ./ reshape(max.([norm(B[:, k]) for k in 1:K], _TINY), 1, :)
        Sm = An' * Bn
        absolute_value && (Sm = abs.(Sm))
        S .*= Sm
    end
    perm = _assignment_max(S)
    per = [S[k, perm[k]] for k in 1:K]
    return (permutation=perm, mean_congruence=mean(per), per_component=per)
end

function _normalize_columns!(A)
    for k in axes(A, 2)
        nrm = norm(@view A[:, k]); nrm > 0 && ((@view A[:, k]) ./= nrm)
    end
    A
end

"""Simulate a sparse Poisson CP tensor with optional covariate effects."""
function simulate_tensor(; n=100, p=20, w=10, K=3, Xcov=nothing,
                         gamma_true=nothing, pi0=0.2, alpha_true=3.0,
                         beta_true=2.0, weights=nothing, sparsity=20.0,
                         seed=42, rng=MersenneTwister(seed))
    unsupervised = gamma_true === false
    if unsupervised
        Xcov = nothing; gamma_true = nothing
    else
        Xcov === nothing && (Xcov = hcat(ones(n), rand(rng, Bernoulli(0.5), n), 20 .+ 60 .* rand(rng, n)))
        gamma_true === nothing && (gamma_true = [[0.0, 2 * rand(rng) - 1, 0.6 * rand(rng) - 0.3] for _ in 1:K])
        length(gamma_true) == K || throw(DimensionMismatch("gamma_true must have length K"))
        all(length(g) == size(Xcov, 2) for g in gamma_true) || throw(DimensionMismatch("gamma vectors must match Xcov"))
    end
    U1 = zeros(n, K)
    for k in 1:K
        multiplier = unsupervised ? ones(n) : exp.(Matrix(Xcov) * gamma_true[k])
        active = rand(rng, n) .> pi0
        for i in findall(active)
            U1[i, k] = rand(rng, Gamma(alpha_true, multiplier[i] / beta_true))
        end
    end
    _normalize_columns!(U1)
    U2 = zeros(p, K); block = fld(p, K)
    for k in 1:K
        lo = (k - 1) * block + 1; hi = k == K ? p : k * block
        U2[lo:hi, k] .= 1 .+ rand(rng, Gamma(2, 0.5), hi - lo + 1)
    end
    _normalize_columns!(U2)
    U3 = zeros(w, K); group = max(1, fld(w, K))
    for k in 1:K
        lo = (k - 1) * group + 1; hi = k == K ? w : min(k * group, w)
        lo <= hi && (U3[lo:hi, k] .= rand(rng, Gamma(2, 1), hi - lo + 1))
    end
    _normalize_columns!(U3)
    weights === nothing && (weights = collect(range(2.0, 0.5; length=K)))
    lambda_true = zeros(n, p, w)
    for k in 1:K, m in 1:w, j in 1:p, i in 1:n
        lambda_true[i, j, m] += weights[k] * U1[i, k] * U2[j, k] * U3[m, k]
    end
    lambda_true .*= sqrt(n * p * w) / sparsity
    X = map(lambda -> rand(rng, Poisson(max(lambda, 0))), lambda_true)
    return (X=X, lambda_true=lambda_true, U1_true=U1, U2_true=U2, U3_true=U3,
            weights=Float64.(weights), Xcov=Xcov, gamma_true=gamma_true,
            sparsity_pct=100 * count(iszero, X) / length(X))
end
