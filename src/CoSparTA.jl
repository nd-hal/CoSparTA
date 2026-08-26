module CoSparTA

using DataFrames
using Distributions
using Dates
using LinearAlgebra
using Optim
using Random
using RecipesBase
using Serialization
using SpecialFunctions
using Statistics
using Tables

export CoSparTAFit, FactorPosterior, PriorEstimate, TensorData, MissingMask
export fit, fit_missing, CoSparTA_missing
export ebpm_point_gamma_multiplier_covariates, ebpm_point_gamma_with_uq, ebps_with_uq
export ebpmf_identity_smooth_control_default
export calc_EZ_3d, calc_EZ_3d_cpp, calc_qz_sparse, calc_qz_sparse_cpp
export adjLF, log_for_ebmf, mKL, poisson_to_multinom
export build_tensor, get_loadings, normalize_factors, project_tensor, reconstruct_tensor
export init_cpapr, select_covariates, match_factors, simulate_tensor
export generate_missing_mask, evaluate_missing_prediction
export get_pip, get_significant_patterns, get_posterior_quantile, get_gamma_ci
export plot_time_factors, plot_channel_factors
export dash_data, dash, dash_stipple

include("types.jl")
include("empirical_bayes.jl")
include("core.jl")
include("preprocessing.jl")
include("postprocessing.jl")
include("missing.jl")
include("inference.jl")
include("visualization.jl")
include("dash.jl")

end
