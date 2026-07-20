const _TINY = 1e-12

_logistic(x::Real) = x >= 0 ? inv(1 + exp(-x)) : exp(x) / (1 + exp(x))
_logit(p::Real) = log(p / (1 - p))

function _logaddexp(a::Real, b::Real)
    m = max(a, b)
    m == Inf && return Inf
    m == -Inf && return -Inf
    return m + log(exp(a - m) + exp(b - m))
end

function _nb_logpdf_logexposure(x::Real, a::Real, logb::Real, logexposure::Real)
    logden = _logaddexp(logb, logexposure)
    gamma_difference = loggamma(a + x) - loggamma(a)
    if x != 0 && gamma_difference == 0
        gamma_difference = digamma(a) * x + trigamma(a) * x^2 / 2
    end
    count_term = iszero(x) ? 0.0 : x * (logexposure - logden)
    return gamma_difference - loggamma(x + 1) +
           a * (logb - logden) + count_term
end

function _finite_hessian(f, x::AbstractVector{<:Real})
    n = length(x)
    H = zeros(Float64, n, n)
    h = 1e-4 .* max.(1.0, abs.(x))
    f0 = f(x)
    for i in 1:n
        xp = Float64.(x); xm = Float64.(x)
        xp[i] += h[i]; xm[i] -= h[i]
        H[i, i] = (f(xp) - 2 * f0 + f(xm)) / h[i]^2
        for j in (i + 1):n
            xpp = Float64.(x); xpm = Float64.(x)
            xmp = Float64.(x); xmm = Float64.(x)
            xpp[i] += h[i]; xpp[j] += h[j]
            xpm[i] += h[i]; xpm[j] -= h[j]
            xmp[i] -= h[i]; xmp[j] += h[j]
            xmm[i] -= h[i]; xmm[j] -= h[j]
            H[i, j] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4 * h[i] * h[j])
            H[j, i] = H[i, j]
        end
    end
    return H
end

function _initial_point_gamma(x, s)
    xmax = max(maximum(x), 1.0)
    smax = max(maximum(s), floatmin(Float64))
    log_rate = log(sum(x ./ xmax) + floatmin(Float64)) + log(xmax) -
               log(sum(s ./ smax)) - log(smax)
    rate_mean = exp(clamp(log_rate, log(1e-6), log(1e6)))
    scaled = clamp.(x ./ max.(s, floatmin(Float64)), 0.0, 1e12)
    v = length(scaled) > 1 ? var(scaled; corrected=false) : rate_mean
    shape = clamp(rate_mean^2 / max(v - rate_mean / max(mean(s), _TINY), 1e-3), 0.05, 100.0)
    beta = clamp(shape / max(rate_mean, 1e-3), 1e-3, 1000.0)
    pi0 = clamp(count(iszero, x) / length(x), 0.01, 0.99)
    return [pi0, shape, beta]
end

function _point_gamma_objective(theta, x, s, X)
    pi0 = _logistic(theta[1])
    a = exp(clamp(theta[2], -20, 20))
    b = exp(clamp(theta[3], -20, 20))
    if !(0.001 <= pi0 <= 0.999 && 0.001 <= a <= 1000 && 0.001 <= b <= 1000)
        return 1e12
    end
    gamma = X === nothing ? nothing : @view theta[4:end]
    ll = 0.0
    for i in eachindex(x)
        linear_predictor = X === nothing ? 0.0 :
                           clamp(dot(@view(X[i, :]), gamma), -100, 100)
        lnb = _nb_logpdf_logexposure(x[i], a, log(b),
                                      log(s[i]) + linear_predictor)
        ll += iszero(x[i]) ? log(pi0 + (1 - pi0) * exp(lnb)) : log1p(-pi0) + lnb
    end
    return isfinite(ll) ? -ll : 1e12
end

function _optimize_pg(x, s, X, g_init; maxiter=500)
    p0 = g_init === nothing ? _initial_point_gamma(x, s) : Float64.(g_init[1:3])
    gamma0 = X === nothing ? Float64[] :
             (g_init === nothing ? zeros(size(X, 2)) : Float64.(g_init[4:end]))
    theta0 = vcat(_logit(clamp(p0[1], 0.001, 0.999)), log(max(p0[2], 0.001)),
                  log(max(p0[3], 0.001)), gamma0)
    objective = t -> _point_gamma_objective(t, x, s, X)
    result = try
        optimize(objective, theta0, LBFGS(), Optim.Options(iterations=maxiter, show_trace=false))
    catch
        optimize(objective, theta0, NelderMead(), Optim.Options(iterations=maxiter, show_trace=false))
    end
    theta = collect(Optim.minimizer(result))
    value = all(isfinite, theta) ? objective(theta) : Inf
    if !isfinite(value) || value >= 1e11
        fallback = optimize(objective, theta0, NelderMead(),
                            Optim.Options(iterations=maxiter, show_trace=false))
        fallback_theta = collect(Optim.minimizer(fallback))
        fallback_value = all(isfinite, fallback_theta) ? objective(fallback_theta) : Inf
        if fallback_value < value
            theta = fallback_theta
            value = fallback_value
        end
    end
    return theta, value, objective
end

function _point_gamma_result(x, s, X; g_init=nothing, maxiter=500, compute_hessian=false)
    xv = Float64.(x)
    sv = length(s) == 1 ? fill(Float64(first(s)), length(xv)) : Float64.(s)
    length(sv) == length(xv) || throw(DimensionMismatch("s must be scalar or match x"))
    all(isfinite, xv) || throw(ArgumentError("x must contain only finite values"))
    all(isfinite, sv) || throw(ArgumentError("s must contain only finite values"))
    any(xv .< 0) && throw(ArgumentError("x must be non-negative"))
    any(sv .<= 0) && throw(ArgumentError("s must be strictly positive"))
    Xm = X === nothing ? nothing : Matrix{Float64}(X)
    Xm !== nothing && size(Xm, 1) != length(xv) && throw(DimensionMismatch("X must have one row per observation"))

    theta, minimum_value, objective = _optimize_pg(xv, sv, Xm, g_init; maxiter=maxiter)
    all(isfinite, theta) || throw(ArgumentError("point-gamma optimization returned non-finite parameters"))
    pi0 = clamp(_logistic(theta[1]), 0.001, 0.999)
    a = exp(clamp(theta[2], log(0.001), log(1000.0)))
    b = exp(clamp(theta[3], log(0.001), log(1000.0)))
    gamma = Xm === nothing ? nothing : clamp.(collect(theta[4:end]), -100.0, 100.0)
    theta_used = vcat(_logit(pi0), log(a), log(b),
                      gamma === nothing ? Float64[] : gamma)
    minimum_value = objective(theta_used)
    multiplier = Xm === nothing ? ones(length(xv)) : exp.(clamp.(Xm * gamma, -100, 100))
    beta_eff = b ./ multiplier
    pi_hat = zeros(length(xv))
    for i in eachindex(xv)
        if iszero(xv[i])
            nb0 = exp(_nb_logpdf_logexposure(0.0, a, log(b),
                                              log(sv[i]) + log(multiplier[i])))
            pi_hat[i] = pi0 / (pi0 + (1 - pi0) * nb0)
        end
    end
    slab_mean = (a .+ xv) ./ (beta_eff .+ sv)
    post_mean = (1 .- pi_hat) .* slab_mean
    post_log = fill(-Inf, length(xv))
    positive = xv .> 0
    post_log[positive] .= (1 .- pi_hat[positive]) .* (digamma.(a .+ xv[positive]) .- log.(beta_eff[positive] .+ sv[positive]))
    post_var = (1 .- pi_hat) .* (a .+ xv) ./ (beta_eff .+ sv).^2 .+
               pi_hat .* (1 .- pi_hat) .* slab_mean.^2
    H = compute_hessian ? _finite_hessian(objective, theta_used) : nothing
    prior = PriorEstimate(pi0, a, inv(b), gamma,
                          Xm === nothing ? :point_gamma : :covariate_dependent, H)
    return EBResult(prior, post_mean, post_log, post_var, 1 .- pi_hat,
                    a .+ xv, beta_eff .+ sv, fill(NaN, length(xv)),
                    -minimum_value, :point_gamma)
end

"""Fit a point-gamma empirical-Bayes Poisson model and return uncertainty summaries."""
function ebpm_point_gamma_with_uq(x, s=1; g_init=nothing, maxiter=500)
    sv = s isa Number ? [Float64(s)] : s
    _point_gamma_result(x, sv, nothing; g_init=g_init, maxiter=maxiter)
end

"""Fit the covariate-dependent point-gamma model used for mode 1."""
function ebpm_point_gamma_multiplier_covariates(x, s, X; g_init=nothing, maxiter=500)
    sv = s isa Number ? [Float64(s)] : s
    _point_gamma_result(x, sv, X; g_init=g_init, maxiter=maxiter, compute_hessian=true)
end

"""
Native Julia replacement for the R `smashrgen::ebps` backend.

It fits a Poisson log-rate with a second-difference Gaussian penalty and uses a
Laplace approximation for posterior log-variance. `level_precision` makes the
otherwise rank-deficient second-difference prior proper, and `max_var_log`
bounds the log-scale Laplace variance in weakly identified components.
"""
function ebps_with_uq(x, s=1; smoothness=10.0, level_precision=1e-3,
                      max_var_log=4.0, maxiter=50, tol=1e-7)
    xv = Float64.(x)
    sv = s isa Number ? fill(Float64(s), length(xv)) : Float64.(s)
    length(sv) == length(xv) || throw(DimensionMismatch("s must be scalar or match x"))
    all(isfinite, xv) || throw(ArgumentError("x must contain only finite values"))
    all(isfinite, sv) || throw(ArgumentError("s must contain only finite values"))
    any(xv .< 0) && throw(ArgumentError("x must be non-negative"))
    any(sv .<= 0) && throw(ArgumentError("s must be strictly positive"))
    smoothness >= 0 || throw(ArgumentError("smoothness must be non-negative"))
    level_precision > 0 || throw(ArgumentError("level_precision must be positive"))
    max_var_log > 0 || throw(ArgumentError("max_var_log must be positive"))
    n = length(xv)
    D = if n >= 3
        M = zeros(n - 2, n)
        for i in 1:(n - 2)
            M[i, i] = 1; M[i, i + 1] = -2; M[i, i + 2] = 1
        end
        M
    else
        zeros(0, n)
    end
    penalty = smoothness .* (D' * D) +
              level_precision .* Matrix{Float64}(I, n, n)
    safe_s = max.(sv, floatmin(Float64))
    log_s = log.(safe_s)
    baseline = log(sum(xv) + 0.5) - log(max(sum(sv), floatmin(Float64)))
    eta = clamp.(log.(xv .+ 0.5) .- log_s, -300.0, 300.0)

    function state(current_eta)
        log_mu = clamp.(log_s .+ current_eta, -700.0, 700.0)
        mu = exp.(log_mu)
        centered = current_eta .- baseline
        objective = sum(mu .- xv .* current_eta .+ loggamma.(xv .+ 1)) +
                    smoothness / 2 * sum(abs2, D * current_eta) +
                    level_precision / 2 * sum(abs2, centered)
        return mu, objective
    end

    for _ in 1:maxiter
        mu, objective = state(eta)
        H = Matrix(Diagonal(mu) + penalty)
        grad = mu .- xv .+ smoothness .* (D' * (D * eta)) .+
               level_precision .* (eta .- baseline)
        step = cholesky(Symmetric(H)) \ grad
        max_step = maximum(abs, step)
        max_step > 10 && (step .*= 10 / max_step)
        eta_new = eta
        step_size = 1.0
        for _ in 1:30
            candidate = clamp.(eta .- step_size .* step, -300.0, 300.0)
            _, candidate_objective = state(candidate)
            if isfinite(candidate_objective) && candidate_objective <= objective
                eta_new = candidate
                break
            end
            step_size /= 2
        end
        maximum(abs.(eta_new .- eta)) < tol && (eta = eta_new; break)
        eta = eta_new
    end
    mu, objective = state(eta)
    H = Matrix(Diagonal(mu) + penalty)
    Hinv = inv(cholesky(Symmetric(H)))
    v_log = clamp.(diag(Hinv), 0.0, max_var_log)
    post_mean = exp.(clamp.(eta .+ v_log ./ 2, -700.0, 700.0))
    post_var = similar(post_mean)
    for i in eachindex(post_var)
        log_variance = 2 * eta[i] + v_log[i] +
                       log(max(expm1(v_log[i]), floatmin(Float64)))
        post_var[i] = exp(clamp(log_variance, -700.0, 700.0))
    end
    prior = PriorEstimate(0.0, NaN, NaN, nothing, :smooth_lognormal, nothing)
    EBResult(prior, post_mean, eta, post_var, fill(NaN, n), fill(NaN, n),
             fill(NaN, n), v_log, -objective, :smooth_lognormal)
end

"""Default options corresponding to the R package's smooth-factor control."""
ebpmf_identity_smooth_control_default() =
    (method=:second_difference, smoothness=10.0, level_precision=1e-3,
     max_var_log=4.0, maxiter=50, tol=1e-7)
