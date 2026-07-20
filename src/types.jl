"""Estimated hyperparameters for one empirical-Bayes component."""
struct PriorEstimate
    pi0::Float64
    shape::Float64
    scale::Float64
    gamma::Union{Nothing,Vector{Float64}}
    kind::Symbol
    hessian::Union{Nothing,Matrix{Float64}}
end

"""Posterior summaries for one tensor mode."""
struct FactorPosterior
    mean::Matrix{Float64}
    mean_log::Matrix{Float64}
    variance::Matrix{Float64}
    pip::Matrix{Float64}
    shape_post::Matrix{Float64}
    rate_post::Matrix{Float64}
    var_log::Matrix{Float64}
    family::Symbol
end

"""Result returned by [`fit`](@ref) and [`fit_missing`](@ref)."""
struct CoSparTAFit
    elbo::Float64
    objective_trace::Vector{Float64}
    U1::FactorPosterior
    U2::FactorPosterior
    U3::FactorPosterior
    U1_normed::FactorPosterior
    U2_normed::FactorPosterior
    U3_normed::FactorPosterior
    weights::Vector{Float64}
    weights_raw::Vector{Float64}
    order::Vector{Int}
    priors1::Vector{Union{Nothing,PriorEstimate}}
    priors2::Vector{Union{Nothing,PriorEstimate}}
    priors3::Vector{Union{Nothing,PriorEstimate}}
    library_size::Vector{Float64}
    observed_mask::Union{Nothing,BitArray{3}}
    runtime_seconds::Float64
    converged::Bool
    iterations::Int
end

"""Output of [`build_tensor`](@ref)."""
struct TensorData{T,L1,L2,L3}
    tensor::Array{T,3}
    dim1_labels::Vector{L1}
    dim2_labels::Vector{L2}
    dim3_labels::Vector{L3}
end

Base.getproperty(x::TensorData, s::Symbol) = s === :X ? getfield(x, :tensor) : getfield(x, s)

"""A held-out tensor and its bookkeeping, returned by [`generate_missing_mask`](@ref)."""
struct MissingMask
    X_obs::Array{Union{Missing,Float64},3}
    obs_mask::BitArray{3}
    obs_indices::Vector{CartesianIndex{3}}
    missing_nonzero_indices::Vector{CartesianIndex{3}}
    true_values::Vector{Float64}
end

Base.getproperty(x::MissingMask, s::Symbol) =
    s === :n_missing ? length(getfield(x, :missing_nonzero_indices)) : getfield(x, s)

struct EBResult
    prior::PriorEstimate
    mean::Vector{Float64}
    mean_log::Vector{Float64}
    variance::Vector{Float64}
    pip::Vector{Float64}
    shape_post::Vector{Float64}
    rate_post::Vector{Float64}
    var_log::Vector{Float64}
    log_likelihood::Float64
    family::Symbol
end

_nanmatrix(n::Integer, k::Integer) = fill(NaN, n, k)

function _posterior(mean, mean_log, variance, pip, shape_post, rate_post, var_log, family)
    FactorPosterior(Matrix{Float64}(mean), Matrix{Float64}(mean_log),
                    Matrix{Float64}(variance), Matrix{Float64}(pip),
                    Matrix{Float64}(shape_post), Matrix{Float64}(rate_post),
                    Matrix{Float64}(var_log), family)
end
