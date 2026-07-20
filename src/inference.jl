function _factor_posterior(result::CoSparTAFit, mode, normalized)
    key = uppercase(String(mode))
    key == "L" && return normalized ? result.U1_normed : result.U1
    key == "F" && return normalized ? result.U2_normed : result.U2
    key == "W" && return normalized ? result.U3_normed : result.U3
    throw(ArgumentError("mode must be L, F, or W"))
end

"""Return posterior inclusion probabilities, optionally thresholded."""
function get_pip(result::CoSparTAFit; mode=:L, threshold=nothing, normalized=true)
    P = _factor_posterior(result, mode, normalized).pip
    all(isnan, P) && return nothing
    if threshold !== nothing
        0 <= threshold <= 1 || throw(ArgumentError("threshold must be between 0 and 1"))
        return P .> threshold
    end
    return P
end

get_pip(result::CoSparTAFit, mode; kwargs...) = get_pip(result; mode=mode, kwargs...)

function _lfdr_discovery(values, alpha)
    valid = findall(isfinite, values)
    isempty(valid) && return Int[]
    order = valid[sortperm(values[valid])]
    cumulative = cumsum(values[order]) ./ (1:length(order))
    accepted = findall(<=(alpha), cumulative)
    isempty(accepted) ? Int[] : sort(order[1:last(accepted)])
end

"""Discover active time points and channels using Bayesian local-FDR control."""
function get_significant_patterns(result::CoSparTAFit; alpha=0.05, mode=:both, normalized=true)
    0 < alpha < 1 || throw(ArgumentError("alpha must lie strictly between zero and one"))
    key = lowercase(String(mode)); key in ("f", "w", "both") || throw(ArgumentError("mode must be F, W, or both"))
    K = size(result.U1.mean, 2)
    Fpip = get_pip(result; mode=:F, normalized=normalized)
    Wpip = get_pip(result; mode=:W, normalized=normalized)
    return [begin
        times = key in ("f", "both") && Fpip !== nothing ? _lfdr_discovery(1 .- Fpip[:, k], alpha) : nothing
        channels = key in ("w", "both") && Wpip !== nothing ? _lfdr_discovery(1 .- Wpip[:, k], alpha) : nothing
        (factor=k, active_times=times, active_channels=channels,
         n_active_times=times === nothing ? 0 : length(times),
         n_active_channels=channels === nothing ? 0 : length(channels))
    end for k in 1:K]
end

function _quantile_name(probability)
    percentage = round(probability * 100; digits=1)
    "q" * (isinteger(percentage) ? string(Int(percentage)) : string(percentage))
end

"""
Compute exact point-gamma or Laplace log-normal posterior quantiles.

The return value is a dictionary keyed like the R result (`"q2.5"`,
`"q97.5"`) with an additional `"runtime_secs"` entry.
"""
function get_posterior_quantile(result::CoSparTAFit; probs=(0.025, 0.975), mode=:L,
                                normalized=true, verbose=false)
    started = time(); P = _factor_posterior(result, mode, normalized)
    all(p -> 0 <= p <= 1, probs) || throw(ArgumentError("probs must lie between zero and one"))
    output = Dict{String,Any}()
    if P.family == :point_gamma
        all(isnan, P.shape_post) && throw(ArgumentError("posterior gamma parameters are unavailable"))
        for tau in probs
            Q = zeros(size(P.mean))
            for idx in eachindex(Q)
                isnan(P.pip[idx]) && (Q[idx] = NaN; continue)
                spike = 1 - P.pip[idx]
                if tau > spike && spike < 1
                    adjusted = clamp((tau - spike) / (1 - spike), 0, 1)
                    Q[idx] = quantile(Gamma(P.shape_post[idx], inv(P.rate_post[idx])), adjusted)
                end
            end
            output[_quantile_name(tau)] = Q
        end
    elseif P.family == :smooth_lognormal
        V = copy(P.var_log)
        missing_var = .!isfinite.(V)
        V[missing_var] .= max.(P.variance[missing_var], 0) ./ (P.mean[missing_var].^2 .+ 1e-10)
        for tau in probs
            output[_quantile_name(tau)] = exp.(P.mean_log .+ sqrt.(max.(V, 0)) .* quantile(Normal(), tau))
        end
    else
        throw(ArgumentError("unknown posterior family $(P.family)"))
    end
    output["runtime_secs"] = time() - started
    verbose && @info "Posterior quantiles computed" runtime_seconds=output["runtime_secs"]
    return output
end

function _gamma_delta(prior::PriorEstimate, level)
    gamma = prior.gamma
    gamma === nothing && return nothing
    H = prior.hessian
    H === nothing && throw(ArgumentError("Hessian unavailable for covariate-dependent component"))
    V = pinv(Symmetric(H)); indices = 4:(3 + length(gamma))
    se = sqrt.(max.(diag(V[indices, indices]), 0))
    z = quantile(Normal(), (1 + level) / 2)
    pvalue = 2 .* ccdf.(Normal(), abs.(gamma ./ max.(se, _TINY)))
    return (estimate=gamma, se=se, lower=gamma .- z .* se, upper=gamma .+ z .* se,
            pvalue=pvalue, method=:delta)
end

"""Confidence intervals for component-specific covariate effects."""
function get_gamma_ci(result::CoSparTAFit; method=:delta, level=0.95, B=200,
                      X=nothing, K=nothing, Xcov=nothing, normalized=true,
                      init_fn=nothing, verbose=true, rng=Random.default_rng(), fit_kwargs...)
    0 < level < 1 || throw(ArgumentError("level must lie strictly between zero and one"))
    priors = normalized ? result.priors1[result.order] : result.priors1
    method = Symbol(method)
    method == :delta && return [p === nothing || p.kind != :covariate_dependent ? nothing : _gamma_delta(p, level) for p in priors]
    method == :bootstrap || throw(ArgumentError("method must be :delta or :bootstrap"))
    (X === nothing || K === nothing || Xcov === nothing) &&
        throw(ArgumentError("X, K, and Xcov are required for bootstrap intervals"))
    lambda = reconstruct_tensor(result)
    draws = [p === nothing || p.gamma === nothing ? nothing : fill(NaN, B, length(p.gamma)) for p in priors]
    reference = normalize_factors(result)
    base_options = merge((verbose=false,), (; fit_kwargs...))
    for b in 1:B
        Xstar = map(rate -> rand(rng, Poisson(max(rate, 0))), lambda)
        options = base_options
        if init_fn !== nothing
            initialization = try init_fn(Xstar, K) catch; nothing end
            initialization !== nothing && (options = merge(options, (init=initialization,)))
        end
        star = try CoSparTA.fit(Xstar, K; Xcov=Xcov, options...) catch exception
            verbose && @warn "Bootstrap fit failed" replicate=b exception=exception
            nothing
        end
        star === nothing && continue
        nf = normalize_factors(star)
        permutation = match_factors((reference.El, reference.Ef, reference.Ew),
                                    (nf.El, nf.Ef, nf.Ew)).permutation
        star_priors = star.priors1[star.order]
        for k in eachindex(draws)
            draws[k] === nothing && continue
            candidate = star_priors[permutation[k]]
            (candidate === nothing || candidate.gamma === nothing) && continue
            length(candidate.gamma) == size(draws[k], 2) && (draws[k][b, :] .= candidate.gamma)
        end
        verbose && b % 10 == 0 && @info "Bootstrap progress" replicate=b total=B
    end
    lo = (1 - level) / 2; hi = (1 + level) / 2
    return map(eachindex(draws)) do k
        draws[k] === nothing && return nothing
        original = priors[k].gamma
        cols = [filter(isfinite, draws[k][:, j]) for j in axes(draws[k], 2)]
        se = [length(c) > 1 ? std(c) : NaN for c in cols]
        lower = [isempty(c) ? NaN : quantile(c, lo) for c in cols]
        upper = [isempty(c) ? NaN : quantile(c, hi) for c in cols]
        pvalue = [isempty(c) ? NaN : min(1.0, 2 * min(mean(c .>= 0), mean(c .<= 0))) for c in cols]
        (estimate=original, se=se, lower=lower, upper=upper, pvalue=pvalue, method=:bootstrap)
    end
end
