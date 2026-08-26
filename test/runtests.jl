using CoSparTA
using DataFrames
using Dates
using LinearAlgebra
using Random
using Serialization
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

@testset "dash_data" begin
    @testset "name-keyed placement (union mode)" begin
        gamma_list = [[0.1, 0.2], [0.3, 0.4]]
        result = dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                            gamma_list=gamma_list, intercept=false,
                            covariate_names=[["age", "income"], ["income", "region"]])
        @test result.gamma_table_mode == :union
        gt = result.gamma_table
        rowA = gt[gt.Rank .== "R1", :]
        rowB = gt[gt.Rank .== "R2", :]
        @test rowA.age[1] == 0.1
        @test rowA.income[1] == 0.2
        @test ismissing(rowA.region[1])
        @test ismissing(rowB.age[1])
        @test rowB.income[1] == 0.3
        @test rowB.region[1] == 0.4
    end

    @testset "heterogeneous overlapping covariate sets (union mode)" begin
        gamma_list = [[0.1, 0.2, 0.3],      # comp1: age, income, region
                      [0.4, 0.5],            # comp2: age, income
                      [0.6, 0.7, 0.8]]       # comp3: income, region, tenure
        covariate_names = [["age", "income", "region"],
                            ["age", "income"],
                            ["income", "region", "tenure"]]
        result = dash_data(Ef=rand(4, 3), Ew=rand(3, 3), lambda=[1.0, 1.0, 1.0],
                            gamma_list=gamma_list, intercept=false,
                            covariate_names=covariate_names)
        @test result.gamma_table_mode == :union
        gt = result.gamma_table
        @test Set(names(gt)) == Set(["Rank", "Note", "age", "income", "region", "tenure"])
        @test length(names(gt)) == 6  # income and region each appear once, not duplicated

        row1 = gt[gt.Rank .== "R1", :]
        row2 = gt[gt.Rank .== "R2", :]
        row3 = gt[gt.Rank .== "R3", :]

        @test row1.age[1] == 0.1
        @test row1.income[1] == 0.2
        @test row1.region[1] == 0.3
        @test ismissing(row1.tenure[1])

        @test row2.age[1] == 0.4
        @test row2.income[1] == 0.5  # comp2's 2nd value must land in "income", not "2nd column"
        @test ismissing(row2.region[1])
        @test ismissing(row2.tenure[1])

        @test ismissing(row3.age[1])
        @test row3.income[1] == 0.6  # comp3's 1st value must land in "income", not "1st column"
        @test row3.region[1] == 0.7
        @test row3.tenure[1] == 0.8
    end

    @testset "disjoint name sets (pairs mode)" begin
        gamma_list = [[0.1, 0.2], [0.3, 0.4]]
        result = dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                            gamma_list=gamma_list, intercept=false,
                            covariate_names=[["age", "tenure"], ["income", "region"]])
        @test result.gamma_table_mode == :pairs
        @test "Covariates" in names(result.gamma_table)
    end

    @testset "unsupervised component note" begin
        gamma_list = [[0.1, 0.2], nothing]
        result = dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                            gamma_list=gamma_list, intercept=false,
                            covariate_names=[["age", "tenure"], nothing])
        gt = result.gamma_table
        row2 = gt[gt.Rank .== "R2", :]
        @test row2.Note[1] == "no covariates for this component"
    end

    @testset "intercept column presence" begin
        with_int = dash_data(Ef=rand(4, 1), Ew=rand(3, 1), lambda=[1.0],
                              gamma_list=[[1.5, 0.1, 0.2]], intercept=true,
                              covariate_names=[["age", "income"]])
        @test "Intercept" in names(with_int.gamma_table)
        @test with_int.covariate_names[1] == ["age", "income"]

        no_int = dash_data(Ef=rand(4, 1), Ew=rand(3, 1), lambda=[1.0],
                            gamma_list=[[0.1, 0.2]], intercept=false,
                            covariate_names=[["age", "income"]])
        @test !("Intercept" in names(no_int.gamma_table))
    end

    @testset "intercept auto-detection" begin
        gamma_list = [[0.1, 0.2, 0.3], nothing]
        result = dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                            gamma_list=gamma_list, covariate_names=["a", "b", "c"])
        @test result.intercept == false
    end

    @testset "vector-form names applied to all supervised components" begin
        gamma_list = [[0.1, 0.2], [0.3, 0.4]]
        result = dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                            gamma_list=gamma_list, intercept=false,
                            covariate_names=["age", "income"])
        @test result.covariate_names[1] == ["age", "income"]
        @test result.covariate_names[2] == ["age", "income"]
        @test result.gamma_table_mode == :union
        gt = result.gamma_table
        @test Set(names(gt)) == Set(["Rank", "Note", "age", "income"])
    end

    @testset "error cases" begin
        @test_throws ErrorException dash_data(Ef=rand(4, 2))
        @test_throws ErrorException dash_data(Ew=rand(3, 2), lambda=[1.0, 1.0])

        gamma_list_uneven = [[0.1, 0.2], [0.3, 0.4, 0.5]]
        @test_throws ErrorException dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                                               gamma_list=gamma_list_uneven, intercept=false,
                                               covariate_names=["a", "b"])

        @test_throws ErrorException dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                                               gamma_list=[[0.1, 0.2], [0.3, 0.4]], intercept=false,
                                               covariate_names=[["a", "b"]])

        @test_throws ErrorException dash_data(Ef=rand(4, 2), Ew=rand(3, 2), lambda=[1.0, 1.0],
                                               gamma_list=[[0.1, 0.2], [0.3, 0.4]], intercept=false,
                                               covariate_names=[["a", "b"], ["c"]])
    end

    @testset "fit path" begin
        sim = simulate_tensor(n=30, p=8, w=6, K=2)
        model = CoSparTA.fit(sim.X, 2; Xcov=sim.Xcov, verbose=false)
        result = dash_data(fit=model)
        @test result.gamma_table isa DataFrame
    end
end

@testset "dash launcher (no server)" begin
    import Pluto

    sim = simulate_tensor(n=30, p=8, w=6, K=2)
    model = CoSparTA.fit(sim.X, 2; Xcov=sim.Xcov, verbose=false)
    pl = dash_data(fit=model)

    dir = mktempdir()
    nb = CoSparTA._dash_prepare(pl, Pluto; dir=dir)

    @test isfile(nb)
    nb_text = read(nb, String)
    @test startswith(nb_text, "### A Pluto.jl notebook ###")
    @test !occursin("{{", nb_text)

    @test try
        open(nb) do io
            Pluto.load_notebook_nobackup(io, nb)
        end
        true
    catch
        try
            # fall back to a path-based signature if the io form is unavailable
            Pluto.load_notebook_nobackup(nb)
            true
        catch
            Pluto.load_notebook(nb)
            true
        end
    end

    payload_path = joinpath(dir, "cosparta_dash_payload.jls")
    @test isfile(payload_path)
    roundtripped = deserialize(payload_path)
    @test roundtripped.K == 2

    @test :dash in names(CoSparTA)
end

@testset "dash_stipple (no server)" begin
    sim = simulate_tensor(n=30, p=8, w=6, K=2)
    model = CoSparTA.fit(sim.X, 2; Xcov=sim.Xcov, verbose=false)
    pl = dash_data(fit=model)

    dir = mktempdir()
    app = CoSparTA._dash_stipple_prepare(pl; dir=dir)

    @test isfile(app)
    app_text = read(app, String)
    @test !occursin("{{PAYLOAD_PATH}}", app_text)
    @test Meta.parseall(app_text) isa Expr

    payload_path = joinpath(dir, "cosparta_dash_payload.jls")
    @test isfile(payload_path)
    roundtripped = deserialize(payload_path)
    @test roundtripped.K == 2

    @test :dash_stipple in names(CoSparTA)
end
