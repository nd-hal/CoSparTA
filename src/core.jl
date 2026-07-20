"""Rescale two factor matrices while preserving their product."""
function adjLF(L::AbstractMatrix, F::AbstractMatrix)
    gamma_l = vec(sum(L; dims=1)); gamma_f = vec(sum(F; dims=1))
    scale = sqrt.(gamma_l .* gamma_f)
    return (L_init=Matrix(L) .* reshape(scale ./ max.(gamma_l, _TINY), 1, :),
            F_init=Matrix(F) .* reshape(scale ./ max.(gamma_f, _TINY), 1, :))
end

"""Library-size-normalized log transform used before empirical-Bayes factorization."""
function log_for_ebmf(Y::AbstractMatrix)
    rs = vec(sum(Y; dims=2))
    med = median(rs)
    return log1p.(med / 0.5 .* Y ./ max.(rs, _TINY))
end

"""Mean generalized Kullback-Leibler divergence."""
function mKL(A, B)
    size(A) == size(B) || throw(DimensionMismatch("A and B must have the same size"))
    vals = Float64[]
    for (a, b) in zip(A, B)
        (ismissing(a) || ismissing(b)) && continue
        av = Float64(a); bv = max(Float64(b), _TINY)
        push!(vals, iszero(av) ? bv : av * log(av / bv) - av + bv)
    end
    return mean(vals)
end

"""Convert a Poisson matrix factorization to multinomial parameterization."""
function poisson_to_multinom(F::AbstractMatrix, L::AbstractMatrix)
    fsums = vec(sum(F; dims=1))
    L2 = Matrix(L) .* reshape(fsums, 1, :)
    s = vec(sum(L2; dims=2))
    L2 ./= max.(s, _TINY)
    F2 = Matrix(F) ./ reshape(max.(fsums, _TINY), 1, :)
    return (FF=F2, L=L2, s=s)
end

"""Sparse responsibility softmax at nonzero tensor coordinates."""
function calc_qz_sparse(I, J, M, Elogl, Elogf, Elogw)
    nnz = length(I); K = size(Elogl, 2)
    length(J) == nnz == length(M) || throw(DimensionMismatch("coordinate lengths differ"))
    alpha = Matrix{Float64}(undef, nnz, K)
    for z in 1:nnz
        rowmax = -Inf
        for k in 1:K
            value = Elogl[I[z], k] + Elogf[J[z], k] + Elogw[M[z], k]
            isnan(value) && (value = -Inf)
            alpha[z, k] = value
            rowmax = max(rowmax, value)
        end
        if rowmax == -Inf
            alpha[z, :] .= inv(K)
            continue
        elseif rowmax == Inf
            infinite = count(==(Inf), @view alpha[z, :])
            for k in 1:K
                alpha[z, k] = alpha[z, k] == Inf ? inv(infinite) : 0.0
            end
            continue
        end
        denom = 0.0
        for k in 1:K
            alpha[z, k] = exp(alpha[z, k] - rowmax)
            denom += alpha[z, k]
        end
        if isfinite(denom) && denom > 0
            alpha[z, :] ./= denom
        else
            alpha[z, :] .= inv(K)
        end
    end
    return alpha
end

calc_qz_sparse_cpp(args...) = calc_qz_sparse(args...)

"""Aggregate weighted sparse counts along all three modes."""
function calc_EZ_3d(I, J, M, values, probabilities, n::Integer, p::Integer, w::Integer)
    rs = zeros(Float64, n); cs = zeros(Float64, p); zs = zeros(Float64, w)
    for z in eachindex(values)
        contribution = values[z] * probabilities[z]
        rs[I[z]] += contribution; cs[J[z]] += contribution; zs[M[z]] += contribution
    end
    return (rs=rs, cs=cs, zs=zs)
end

calc_EZ_3d_cpp(args...) = calc_EZ_3d(args...)

function _coordinates(X, mask)
    I = Int[]; J = Int[]; M = Int[]; values = Float64[]
    for idx in CartesianIndices(X)
        if mask[idx] && X[idx] > 0
            push!(I, idx[1]); push!(J, idx[2]); push!(M, idx[3]); push!(values, X[idx])
        end
    end
    return I, J, M, values
end

function _clean_tensor(Xin, obs_mask)
    ndims(Xin) == 3 || throw(DimensionMismatch("X must be a three-dimensional array"))
    mask = obs_mask === nothing ? trues(size(Xin)) : BitArray(obs_mask)
    size(mask) == size(Xin) || throw(DimensionMismatch("obs_mask must match X"))
    X = zeros(Float64, size(Xin))
    for idx in eachindex(Xin)
        value = Xin[idx]
        if mask[idx] && !ismissing(value)
            v = Float64(value)
            isfinite(v) || throw(ArgumentError("observed values must be finite"))
            v >= 0 || throw(ArgumentError("X must be non-negative"))
            X[idx] = v
        else
            mask[idx] = false
        end
    end
    return X, mask
end

function _column_change(A, B)
    out = 0.0
    for k in axes(A, 2)
        na = norm(@view A[:, k]); nb = norm(@view B[:, k])
        na = na > 0 ? na : 1.0; nb = nb > 0 ? nb : 1.0
        out = max(out, maximum(abs.(@view(A[:, k])) ./ na .- abs.(@view(B[:, k])) ./ nb))
    end
    return out
end

function _reconstruction_change(curr, prev)
    L, F, W = curr; L0, F0, W0 = prev
    denom2 = max(sum((L0' * L0) .* (F0' * F0) .* (W0' * W0)), 1e-20)
    changes = zeros(size(L, 2))
    for k in axes(L, 2)
        nc = sum(abs2, @view L[:, k]) * sum(abs2, @view F[:, k]) * sum(abs2, @view W[:, k])
        np = sum(abs2, @view L0[:, k]) * sum(abs2, @view F0[:, k]) * sum(abs2, @view W0[:, k])
        np < 1e-20 && continue
        dp = dot(@view(L[:, k]), @view(L0[:, k])) * dot(@view(F[:, k]), @view(F0[:, k])) * dot(@view(W[:, k]), @view(W0[:, k]))
        changes[k] = sqrt(max(nc - 2dp + np, 0)) / sqrt(denom2)
    end
    return maximum(changes)
end

function _observed_l_scales(F, W, mask, k)
    n, p, w = size(mask)
    ls = zeros(n)
    for idx in CartesianIndices(mask)
        mask[idx] || continue
        i, j, m = Tuple(idx)
        ls[i] += F[j, k] * W[m, k]
    end
    return max.(ls, 1e-10)
end

function _observed_f_scales(L, W, mask, k)
    n, p, w = size(mask)
    fs = zeros(p)
    for idx in CartesianIndices(mask)
        mask[idx] || continue
        i, j, m = Tuple(idx)
        fs[j] += L[i, k] * W[m, k]
    end
    return max.(fs, 1e-10)
end

function _observed_w_scales(L, F, mask, k)
    n, p, w = size(mask)
    ws = zeros(w)
    for idx in CartesianIndices(mask)
        mask[idx] || continue
        i, j, m = Tuple(idx)
        ws[m] += L[i, k] * F[j, k]
    end
    return max.(ws, 1e-10)
end

function _store!(mean, elog, variance, pip, shape, rate, varlog, priors, k, eb)
    all(isfinite, eb.mean) || throw(ArgumentError("empirical-Bayes update produced non-finite posterior means for component $k"))
    all(x -> isfinite(x) || x == -Inf, eb.mean_log) ||
        throw(ArgumentError("empirical-Bayes update produced invalid expected logs for component $k"))
    mean[:, k] .= eb.mean; elog[:, k] .= eb.mean_log
    variance[:, k] .= eb.variance; pip[:, k] .= eb.pip
    shape[:, k] .= eb.shape_post; rate[:, k] .= eb.rate_post; varlog[:, k] .= eb.var_log
    priors[k] = eb.prior
end

function _rescale_component!(mean, elog, variance, rate, priors, k, factor)
    isfinite(factor) && factor > 0 || return
    mean[:, k] .*= factor
    elog[:, k] .+= log(factor)
    variance[:, k] .*= factor^2
    rate[:, k] ./= factor
    prior = priors[k]
    if prior !== nothing
        priors[k] = PriorEstimate(prior.pi0, prior.shape, prior.scale * factor,
                                  prior.gamma, prior.kind, prior.hessian)
    end
end

function _balance_component!(L, ElogL, VL, RL, priorsL,
                             F, ElogF, VF, RF, priorsF,
                             W, ElogW, VW, RW, priorsW, k,
                             fix_L, fix_F, fix_W)
    arrays = ((L, ElogL, VL, RL, priorsL),
              (F, ElogF, VF, RF, priorsF),
              (W, ElogW, VW, RW, priorsW))
    fixed = (fix_L, fix_F, fix_W)
    mutable_modes = findall(!, collect(fixed))
    length(mutable_modes) >= 2 || return
    sums = [sum(@view arrays[m][1][:, k]) for m in mutable_modes]
    all(x -> isfinite(x) && x > 0, sums) || return
    logtarget = mean(log, sums)
    factors = exp.(logtarget .- log.(sums))
    all(x -> isfinite(x) && x > 0, factors) || return
    for (position, mode) in enumerate(mutable_modes)
        _rescale_component!(arrays[mode]..., k, factors[position])
    end
end

function _expand(M, keep, n_original; fillvalue=0.0)
    out = fill(fillvalue, n_original, size(M, 2)); out[keep, :] .= M; out
end

function _norm_factor(P::FactorPosterior, norms, order)
    scale = reshape(max.(norms, _TINY), 1, :)
    FactorPosterior((P.mean ./ scale)[:, order],
                    (P.mean_log .- log.(scale))[:, order],
                    (P.variance ./ scale.^2)[:, order], P.pip[:, order],
                    P.shape_post[:, order], (P.rate_post .* scale)[:, order],
                    P.var_log[:, order], P.family)
end

function _elbo_sparse(I, J, M, values, alpha, L, F, W, mask)
    value = 0.0
    for z in eachindex(values), k in axes(L, 2)
        a = max(alpha[z, k], _TINY)
        rate = max(L[I[z], k] * F[J[z], k] * W[M[z], k], _TINY)
        value += values[z] * a * (log(rate) - log(a))
    end
    total_rate = 0.0
    if all(mask)
        for k in axes(L, 2)
            total_rate += sum(@view(L[:, k])) * sum(@view(F[:, k])) * sum(@view(W[:, k]))
        end
    else
        for idx in CartesianIndices(mask)
            mask[idx] || continue
            i, j, m = Tuple(idx)
            for k in axes(L, 2)
                total_rate += L[i, k] * F[j, k] * W[m, k]
            end
        end
    end
    return value - total_rate - sum(loggamma(v + 1) for v in values)
end

function _resolve_xcov(Xcov, K, n_original, keep)
    Xcov === nothing && return fill(nothing, K)
    xs = Xcov isa AbstractMatrix ? [Xcov for _ in 1:K] : collect(Xcov)
    length(xs) == K || throw(DimensionMismatch("Xcov must be a matrix or a length-K collection"))
    return map(xs) do x
        x === nothing && return nothing
        size(x, 1) == n_original || throw(DimensionMismatch("each Xcov matrix must have size(X,1) rows"))
        out = Matrix{Float64}(x[keep, :])
        all(isfinite, out) || throw(ArgumentError("Xcov must contain only finite values"))
        out
    end
end

function _fit_impl(Xin, K::Integer; Xcov=nothing, lib_size=nothing,
                   init=:random_gamma, maxiter=100, maxiter_init=100, tol=1e-6,
                   compute_elbo_final=false, fix_L=false, fix_F=false, fix_W=false,
                   smooth_F=true, printevery=10, verbose=true,
                   convergence_criteria=:factor_change, n_stable=3,
                   ebpm_fns=nothing, obs_mask=nothing, missing_mode=false,
                   U1_true=nothing, U2_true=nothing, U3_true=nothing,
                   rng=Random.default_rng())
    K > 0 || throw(ArgumentError("K must be positive"))
    start_time = time()
    Xall, maskall = _clean_tensor(Xin, obs_mask)
    original_size = size(Xall)
    user_keep = vec(dropdims(sum(Xall; dims=(2, 3)); dims=(2, 3))) .> 0
    time_keep = vec(dropdims(sum(Xall; dims=(1, 3)); dims=(1, 3))) .> 0
    channel_keep = vec(dropdims(sum(Xall; dims=(1, 2)); dims=(1, 2))) .> 0
    all(user_keep) || verbose && @info "Dropping $(count(!, user_keep)) all-zero mode-1 slices"
    any(user_keep) && any(time_keep) && any(channel_keep) || throw(ArgumentError("X has no nonzero three-way support"))
    X = Xall[user_keep, time_keep, channel_keep]
    mask = maskall[user_keep, time_keep, channel_keep]
    n, p, w = size(X)
    xcovs = _resolve_xcov(Xcov, K, original_size[1], user_keep)
    libs = lib_size === nothing ? ones(n) : begin
        length(lib_size) == original_size[1] || throw(DimensionMismatch("lib_size must have size(X,1) entries"))
        out = Float64.(lib_size[user_keep])
        all(x -> isfinite(x) && x > 0, out) ||
            throw(ArgumentError("lib_size must contain finite, strictly positive values"))
        out
    end

    if init == :random_gamma || init == "random_gamma"
        L = rand(rng, Gamma(100, 0.01), n, K)
        F = rand(rng, Gamma(100, 0.01), p, K)
        W = rand(rng, Gamma(100, 0.01), w, K)
    elseif init isa Tuple || init isa AbstractVector
        length(init) == 3 || throw(ArgumentError("init must contain L, F, and W"))
        L = Matrix{Float64}(init[1][user_keep, :])
        F = Matrix{Float64}(init[2][time_keep, :])
        W = Matrix{Float64}(init[3][channel_keep, :])
        size(L, 2) == K == size(F, 2) == size(W, 2) || throw(DimensionMismatch("initial factors must have K columns"))
        all(A -> all(x -> isfinite(x) && x >= 0, A), (L, F, W)) ||
            throw(ArgumentError("initial factors must be finite and non-negative"))
        all(A -> all(vec(sum(A; dims=1)) .> 0), (L, F, W)) ||
            throw(ArgumentError("every initial factor column must contain a positive value"))
    else
        throw(ArgumentError("init must be :random_gamma or a three-factor collection"))
    end
    ElogL = log.(L .+ 1e-18); ElogF = log.(F .+ 1e-18); ElogW = log.(W .+ 1e-18)
    VL = _nanmatrix(n, K); VF = _nanmatrix(p, K); VW = _nanmatrix(w, K)
    PL = _nanmatrix(n, K); PF = _nanmatrix(p, K); PW = _nanmatrix(w, K)
    SL = _nanmatrix(n, K); SF = _nanmatrix(p, K); SW = _nanmatrix(w, K)
    RL = _nanmatrix(n, K); RF = _nanmatrix(p, K); RW = _nanmatrix(w, K)
    VlogL = _nanmatrix(n, K); VlogF = _nanmatrix(p, K); VlogW = _nanmatrix(w, K)
    priorsL = Union{Nothing,PriorEstimate}[nothing for _ in 1:K]
    priorsF = Union{Nothing,PriorEstimate}[nothing for _ in 1:K]
    priorsW = Union{Nothing,PriorEstimate}[nothing for _ in 1:K]
    I, J, M, values = _coordinates(X, mask)
    alpha = calc_qz_sparse(I, J, M, ElogL, ElogF, ElogW)
    trace = Float64[]; stable = 0; converged = false; iterations = 0
    if ebpm_fns === nothing
        l_fns = nothing
        f_fn = smooth_F ? ebps_with_uq : ebpm_point_gamma_with_uq
        w_fn = ebpm_point_gamma_with_uq
    else
        length(ebpm_fns) == 3 || throw(ArgumentError("ebpm_fns must contain L, F, and W functions"))
        l_spec, f_fn, w_fn = ebpm_fns
        l_fns = l_spec isa Function ? [l_spec for _ in 1:K] : collect(l_spec)
        length(l_fns) == K || throw(ArgumentError("the L-mode ebpm_fns entry must be one function or K functions"))
    end
    familyL = :point_gamma; familyF = smooth_F ? :smooth_lognormal : :point_gamma; familyW = :point_gamma

    verbose && @info "Running CoSparTA" dimensions=(n, p, w) rank=K nonzeros=length(values)
    for iter in 1:maxiter
        iterations = iter
        previous = (copy(L), copy(F), copy(W))
        for k in 1:K
            Ez = calc_EZ_3d(I, J, M, values, @view(alpha[:, k]), n, p, w)
            if !fix_L
                ls = if missing_mode
                    _observed_l_scales(F, W, mask, k) .* libs
                else
                    fill(sum(@view(F[:, k])) * sum(@view(W[:, k])), n) .* libs
                end
                eb = if l_fns === nothing
                    xcovs[k] === nothing ? ebpm_point_gamma_with_uq(Ez.rs, ls) :
                        ebpm_point_gamma_multiplier_covariates(Ez.rs, ls, xcovs[k])
                else
                    xcovs[k] === nothing ? l_fns[k](Ez.rs, ls) : l_fns[k](Ez.rs, ls, xcovs[k])
                end
                _store!(L, ElogL, VL, PL, SL, RL, VlogL, priorsL, k, eb)
                familyL = eb.family
            end
            if !fix_F
                fs = missing_mode ? _observed_f_scales(L, W, mask, k) :
                     fill(sum(@view(L[:, k])) * sum(@view(W[:, k])), p)
                eb = f_fn(Ez.cs, fs)
                _store!(F, ElogF, VF, PF, SF, RF, VlogF, priorsF, k, eb)
                familyF = eb.family
            end
            if !fix_W
                ws = missing_mode ? _observed_w_scales(L, F, mask, k) :
                     fill(sum(@view(L[:, k])) * sum(@view(F[:, k])), w)
                eb = w_fn(Ez.zs, ws)
                _store!(W, ElogW, VW, PW, SW, RW, VlogW, priorsW, k, eb)
                familyW = eb.family
            end
            _balance_component!(L, ElogL, VL, RL, priorsL,
                                F, ElogF, VF, RF, priorsF,
                                W, ElogW, VW, RW, priorsW, k,
                                fix_L, fix_F, fix_W)
        end
        alpha = calc_qz_sparse(I, J, M, ElogL, ElogF, ElogW)
        criterion = Symbol(convergence_criteria)
        objective = if criterion == :factor_change
            max(_column_change(L, previous[1]), _column_change(F, previous[2]), _column_change(W, previous[3]))
        elseif criterion == :recon_change
            _reconstruction_change((L, F, W), previous)
        elseif criterion == :mKLabs
            pred = [sum(L[I[z], k] * F[J[z], k] * W[M[z], k] for k in 1:K) for z in eachindex(values)]
            mKL(values, pred)
        elseif criterion == :ELBO
            _elbo_sparse(I, J, M, values, alpha, L, F, W, mask)
        else
            throw(ArgumentError("unknown convergence_criteria: $convergence_criteria"))
        end
        push!(trace, objective)
        if criterion in (:factor_change, :recon_change)
            stable = objective < tol ? stable + 1 : 0
            converged = stable >= n_stable
        elseif length(trace) > 1
            delta = criterion == :ELBO ? (trace[end] - trace[end - 1]) / length(X) : abs(trace[end] - trace[end - 1])
            converged = delta < tol
        end
        verbose && printevery > 0 && iter % printevery == 0 && @info "CoSparTA iteration" iteration=iter objective=objective
        converged && break
    end

    elbo = convergence_criteria in (:ELBO, "ELBO") ? last(trace) :
           (compute_elbo_final ? _elbo_sparse(I, J, M, values, alpha, L, F, W, mask) : NaN)
    Lp = _posterior(_expand(L, user_keep, original_size[1]), _expand(ElogL, user_keep, original_size[1]; fillvalue=-Inf),
                    _expand(VL, user_keep, original_size[1]; fillvalue=NaN), _expand(PL, user_keep, original_size[1]; fillvalue=NaN),
                    _expand(SL, user_keep, original_size[1]; fillvalue=NaN), _expand(RL, user_keep, original_size[1]; fillvalue=NaN),
                    _expand(VlogL, user_keep, original_size[1]; fillvalue=NaN), familyL)
    Fp = _posterior(_expand(F, time_keep, original_size[2]), _expand(ElogF, time_keep, original_size[2]; fillvalue=-Inf),
                    _expand(VF, time_keep, original_size[2]; fillvalue=NaN), _expand(PF, time_keep, original_size[2]; fillvalue=NaN),
                    _expand(SF, time_keep, original_size[2]; fillvalue=NaN), _expand(RF, time_keep, original_size[2]; fillvalue=NaN),
                    _expand(VlogF, time_keep, original_size[2]; fillvalue=NaN), familyF)
    Wp = _posterior(_expand(W, channel_keep, original_size[3]), _expand(ElogW, channel_keep, original_size[3]; fillvalue=-Inf),
                    _expand(VW, channel_keep, original_size[3]; fillvalue=NaN), _expand(PW, channel_keep, original_size[3]; fillvalue=NaN),
                    _expand(SW, channel_keep, original_size[3]; fillvalue=NaN), _expand(RW, channel_keep, original_size[3]; fillvalue=NaN),
                    _expand(VlogW, channel_keep, original_size[3]; fillvalue=NaN), familyW)
    nl = [norm(filter(isfinite, Lp.mean[:, k])) for k in 1:K]
    nf = [norm(filter(isfinite, Fp.mean[:, k])) for k in 1:K]
    nw = [norm(filter(isfinite, Wp.mean[:, k])) for k in 1:K]
    raw_weights = nl .* nf .* nw; order = sortperm(raw_weights; rev=true)
    Ln = _norm_factor(Lp, nl, order); Fn = _norm_factor(Fp, nf, order); Wn = _norm_factor(Wp, nw, order)
    return CoSparTAFit(elbo, trace, Lp, Fp, Wp, Ln, Fn, Wn, raw_weights[order], raw_weights,
                       order, priorsL, priorsF, priorsW, Float64.(lib_size === nothing ? ones(original_size[1]) : lib_size),
                       missing_mode ? maskall : nothing, time() - start_time, converged, iterations)
end

"""Fit the covariate-aware sparse Poisson CP model."""
fit(X, K::Integer; kwargs...) = _fit_impl(X, K; kwargs...)
fit(X; K::Integer, kwargs...) = fit(X, K; kwargs...)

"""Fit CoSparTA while excluding entries marked missing from every update."""
fit_missing(X, K::Integer; obs_mask=nothing, kwargs...) =
    _fit_impl(X, K; obs_mask=obs_mask, missing_mode=true, kwargs...)
fit_missing(X; K::Integer, kwargs...) = fit_missing(X, K; kwargs...)

CoSparTA_missing(args...; kwargs...) = fit_missing(args...; kwargs...)
