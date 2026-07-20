using CoSparTA
using Dates
using LinearAlgebra
using Random
using Test

@testset "sparse kernels" begin
    I = [1, 1, 2]; J = [1, 2, 1]; M = [1, 1, 2]; values = [2.0, 3.0, 4.0]
    E1 = log.([1.0 2.0; 2.0 1.0])
    E2 = log.([1.0 1.0; 2.0 1.0])
    E3 = log.([1.0 2.0; 1.0 1.0])
    q = calc_qz_sparse(I, J, M, E1, E2, E3)
    @test size(q) == (3, 2)
    @test all(isapprox.(vec(sum(q; dims=2)), 1.0; atol=1e-12))
    ez = calc_EZ_3d(I, J, M, values, q[:, 1], 2, 2, 2)
    @test sum(ez.rs) ≈ sum(values .* q[:, 1])
    @test calc_qz_sparse_cpp(I, J, M, E1, E2, E3) ≈ q

    fallback = calc_qz_sparse([1], [1], [1], fill(-Inf, 1, 2),
                              fill(-Inf, 1, 2), fill(-Inf, 1, 2))
    @test fallback == fill(0.5, 1, 2)
end

@testset "preprocessing" begin
    table = (user=["u1", "u1", "u2"], time=[2, 1, 1], channel=["a", "a", "b"])
    built = build_tensor(table; row=:user, col=:time, slice=:channel, verbose=false)
    @test size(built.X) == (2, 2, 2)
    @test sum(built.X) == 3
    @test built.dim2_labels == [1, 2]

    dated = (user=[1, 1, 2], time=DateTime(2024, 1, 1) .+ Hour.([0, 2, 3]), channel=[:a, :a, :b])
    dated_tensor = build_tensor(dated; row=:user, col=:time, slice=:channel,
                                time_bins=2, verbose=false)
    @test size(dated_tensor.X, 2) == 2
end

@testset "core fit and postprocessing" begin
    rng = MersenneTwister(11)
    X = rand(rng, 0:3, 24, 8, 5)
    model = fit(X, 2; maxiter=3, verbose=false, rng=MersenneTwister(2))
    @test size(get_loadings(model, "U1")) == (24, 2)
    @test size(get_loadings(model, "U2")) == (8, 2)
    @test size(get_loadings(model, "U3")) == (5, 2)
    @test all(get_loadings(model, "weight") .> 0)
    @test issorted(get_loadings(model, "weight"); rev=true)
    @test all(abs.(vec(sqrt.(sum(get_loadings(model, "U1").^2; dims=1))) .- 1) .< 1e-8)

    raw = reconstruct_tensor(model; normalized=false)
    normalized = reconstruct_tensor(model)
    @test size(raw) == size(X)
    @test maximum(abs.(raw .- normalized)) < 1e-8
    @test all(normalized .>= 0)

    projected = project_tensor(X, model)
    @test size(projected) == (24, 2)
    @test project_tensor(X[1, :, :], model) ≈ vec(projected[1, :])

    qL = get_posterior_quantile(model; probs=(0.025, 0.975), mode=:L)
    @test size(qL["q2.5"]) == (24, 2)
    @test all(qL["q2.5"] .<= qL["q97.5"])
    qF = get_posterior_quantile(model; mode=:F)
    @test all(qF["q2.5"] .> 0)

    matched = match_factors((model.U1_normed.mean, model.U2_normed.mean),
                            (model.U1_normed.mean[:, [2, 1]], model.U2_normed.mean[:, [2, 1]]))
    @test matched.permutation == [2, 1]
    @test matched.mean_congruence ≈ 1

    screen = select_covariates(K=2, covariate_data=hcat(collect(1.0:24.0), randn(rng, 24)),
                               El=model.U1_normed.mean, verbose=false)
    @test length(screen.selected) == 2
    @test occursin("TimeFactorPlot", string(typeof(plot_time_factors(model))))
    @test occursin("ChannelFactorPlot", string(typeof(plot_channel_factors(model))))
end

@testset "covariates and missing data" begin
    rng = MersenneTwister(21)
    X = rand(rng, 0:3, 18, 6, 4)
    Xcov = hcat(ones(18), randn(rng, 18))
    model = fit(X, 2; Xcov=Xcov, maxiter=1, verbose=false, rng=MersenneTwister(3))
    @test all(p -> p !== nothing && p.gamma !== nothing && length(p.gamma) == 2, model.priors1)
    @test size(get_pip(model; mode=:L)) == (18, 2)
    @test get_pip(model; mode=:F) === nothing
    ci = get_gamma_ci(model; method=:delta)
    @test length(ci) == 2
    @test length(ci[1].estimate) == 2

    heldout = generate_missing_mask(X; missing_rate=0.1, seed=9, verbose=false)
    missing_model = fit_missing(heldout.X_obs, 2; obs_mask=heldout.obs_mask,
                                maxiter=2, verbose=false, rng=MersenneTwister(4),
                                convergence_criteria=:factor_change)
    metrics = evaluate_missing_prediction(missing_model, heldout; verbose=false)
    @test isfinite(metrics.rmse)
    @test isfinite(metrics.mae)
    @test size(reconstruct_tensor(missing_model)) == size(X)
end

@testset "simulation and warm start" begin
    sim = simulate_tensor(n=20, p=8, w=6, K=2, seed=5)
    @test size(sim.X) == (20, 8, 6)
    @test size(sim.U1_true) == (20, 2)
    initialization = init_cpapr(sim.X, 2; n_iters=2, random_state=5)
    @test map(size, initialization) == ((20, 2), (8, 2), (6, 2))
end

@testset "numerical stability regressions" begin
    zero_smooth = ebps_with_uq(zeros(20), ones(20))
    @test all(isfinite, zero_smooth.mean)
    @test all(isfinite, zero_smooth.variance)
    @test all(isfinite, zero_smooth.var_log)
    tiny_smooth = ebps_with_uq(10.0 .^ range(-80, -12; length=20), fill(1e-80, 20))
    @test all(isfinite, tiny_smooth.mean)
    @test all(isfinite, tiny_smooth.variance)
    @test_throws ArgumentError ebps_with_uq(zeros(3), zeros(3))

    # The R algorithm is Gauss-Seidel: F must use the newly updated L, and W
    # must use both the newly updated L and F. This guards against stale scales.
    seen_scales = Dict{Symbol,Vector{Float64}}()
    function fake_eb(x, value)
        n = length(x)
        prior = PriorEstimate(0.0, 1.0, 1.0, nothing, :point_gamma, nothing)
        CoSparTA.EBResult(prior, fill(value, n), fill(log(value), n), zeros(n),
                         ones(n), ones(n), ones(n), fill(NaN, n), 0.0,
                         :point_gamma)
    end
    l_fn = function (x, s)
        seen_scales[:L] = copy(s)
        fake_eb(x, 2.0)
    end
    f_fn = function (x, s)
        seen_scales[:F] = copy(s)
        fake_eb(x, 3.0)
    end
    w_fn = function (x, s)
        seen_scales[:W] = copy(s)
        fake_eb(x, 4.0)
    end
    fit(ones(Int, 2, 3, 2), 1; init=(ones(2, 1), ones(3, 1), ones(2, 1)),
        maxiter=1, smooth_F=false, ebpm_fns=(l_fn, f_fn, w_fn), verbose=false)
    @test seen_scales[:L] == fill(6.0, 2)
    @test seen_scales[:F] == fill(8.0, 3)
    @test seen_scales[:W] == fill(36.0, 2)

    # This is the documented quick-start configuration. The stale-scale porting
    # bug used to overflow the smooth-factor posterior by iteration 6.
    demo = simulate_tensor(n=100, p=20, w=10, K=3, seed=42)
    demo_model = fit(demo.X, 3; Xcov=demo.Xcov, maxiter=8,
                     convergence_criteria=:factor_change,
                     rng=MersenneTwister(1), verbose=false)
    @test all(isfinite, demo_model.U1.mean)
    @test all(isfinite, demo_model.U2.mean)
    @test all(isfinite, demo_model.U3.mean)
    @test all(isfinite, reconstruct_tensor(demo_model))
end
