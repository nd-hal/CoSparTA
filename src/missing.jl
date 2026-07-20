"""Randomly hold out a proportion of nonzero tensor entries."""
function generate_missing_mask(X::AbstractArray; missing_rate=0.1, seed=nothing,
                               rng=seed === nothing ? Random.default_rng() : MersenneTwister(seed),
                               verbose=true)
    0 <= missing_rate <= 1 || throw(ArgumentError("missing_rate must be between 0 and 1"))
    ndims(X) == 3 || throw(DimensionMismatch("X must be three-dimensional"))
    nonzero = [idx for idx in CartesianIndices(X) if !ismissing(X[idx]) && X[idx] > 0]
    nmissing = round(Int, length(nonzero) * missing_rate)
    selected = nmissing == 0 ? Int[] : randperm(rng, length(nonzero))[1:nmissing]
    missing_indices = nonzero[selected]
    true_values = Float64[X[idx] for idx in missing_indices]
    Xobs = Array{Union{Missing,Float64}}(undef, size(X))
    for idx in CartesianIndices(X)
        Xobs[idx] = ismissing(X[idx]) ? missing : Float64(X[idx])
    end
    mask = trues(size(X))
    for idx in missing_indices
        Xobs[idx] = missing; mask[idx] = false
    end
    observed = [idx for idx in CartesianIndices(mask) if mask[idx]]
    verbose && @info "Generated missing mask" total=length(X) original_nonzeros=length(nonzero) held_out=nmissing
    MissingMask(Xobs, mask, observed, missing_indices, true_values)
end

generate_missing_mask(X, missing_rate::Real; kwargs...) =
    generate_missing_mask(X; missing_rate=missing_rate, kwargs...)

"""Evaluate reconstruction at entries held out by [`generate_missing_mask`](@ref)."""
function evaluate_missing_prediction(result::CoSparTAFit, info::MissingMask; verbose=true)
    L = result.U1.mean; F = result.U2.mean; W = result.U3.mean
    prediction = zeros(length(info.missing_nonzero_indices))
    for (z, idx) in enumerate(info.missing_nonzero_indices)
        i, j, m = Tuple(idx)
        prediction[z] = sum(L[i, k] * F[j, k] * W[m, k] for k in axes(L, 2))
    end
    truth = info.true_values
    rmse = isempty(truth) ? NaN : sqrt(mean((prediction .- truth).^2))
    mae = isempty(truth) ? NaN : mean(abs.(prediction .- truth))
    correlation = length(truth) < 2 || iszero(std(prediction)) || iszero(std(truth)) ? NaN : cor(prediction, truth)
    deviance = 2 * sum(t * log(max(t, 1e-10) / max(y, 1e-10)) - (t - y)
                       for (t, y) in zip(truth, prediction))
    verbose && @info "Missing-data prediction" rmse=rmse mae=mae correlation=correlation deviance=deviance
    return (rmse=rmse, mae=mae, correlation=correlation, deviance=deviance,
            predicted=prediction, true_values=truth)
end
