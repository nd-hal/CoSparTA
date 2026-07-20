using CoSparTA
using Random

sim = simulate_tensor(n=100, p=20, w=10, K=3, seed=42)

model = fit(sim.X, 3;
            Xcov=sim.Xcov,
            maxiter=20,
            convergence_criteria=:factor_change,
            tol=1e-6,
            rng=MersenneTwister(42))

println("component weights: ", round.(get_loadings(model, "weight"); digits=3))
println("U1 size: ", size(get_loadings(model, "U1")))
println("U2 size: ", size(get_loadings(model, "U2")))
println("U3 size: ", size(get_loadings(model, "U3")))

ci = get_posterior_quantile(model; mode=:W)
println("channel CI size: ", size(ci["q2.5"]))

heldout = generate_missing_mask(sim.X; missing_rate=0.10, seed=7)
missing_model = fit_missing(heldout.X_obs, 3;
                            Xcov=sim.Xcov,
                            obs_mask=heldout.obs_mask,
                            maxiter=20,
                            convergence_criteria=:factor_change,
                            verbose=false,
                            rng=MersenneTwister(43))
println(evaluate_missing_prediction(missing_model, heldout))
