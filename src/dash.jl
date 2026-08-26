_dash_gamma_len(g) = g === nothing ? 0 : length(g)

function _dash_resolve_covariate_names(covariate_names, gamma_lens, supervised, K,
                                        intercept_flag, intercept_specified)
    if covariate_names === nothing
        return (fill(nothing, K), intercept_flag)
    end

    is_vector_form = covariate_names isa AbstractVector{<:AbstractString}
    is_list_form = !is_vector_form && covariate_names isa AbstractVector &&
        all(x -> x === nothing || x isa AbstractVector{<:AbstractString}, covariate_names)

    if !is_vector_form && !is_list_form
        error("covariate_names must be a vector of strings or a length-K list of per-component names.")
    end

    names = Vector{Union{Nothing,Vector{String}}}(fill(nothing, K))

    if is_list_form
        length(covariate_names) == K ||
            error("covariate_names supplied as a list must have length K=$K; got $(length(covariate_names)).")

        if !intercept_specified
            named_k = [k for k in 1:K if supervised[k] && covariate_names[k] !== nothing]
            if !isempty(named_k)
                matches_no_intercept = all(k -> length(covariate_names[k]) == gamma_lens[k] - 1, named_k)
                matches_with_intercept = all(k -> length(covariate_names[k]) == gamma_lens[k], named_k)
                if !matches_no_intercept && matches_with_intercept
                    intercept_flag = false
                    @warn "covariate_names lengths match gamma lengths exactly, implying no intercept; proceeding with intercept=false. Pass intercept=false explicitly to silence this warning."
                end
            end
        end

        for k in 1:K
            if supervised[k] && covariate_names[k] !== nothing
                expected = intercept_flag ? gamma_lens[k] - 1 : gamma_lens[k]
                length(covariate_names[k]) == expected ||
                    error("covariate_names for component $k must have length $expected; got $(length(covariate_names[k])).")
                names[k] = collect(String, covariate_names[k])
            end
        end

        return (names, intercept_flag)
    end

    supervised_lens = unique(gamma_lens[supervised])
    if length(supervised_lens) > 1
        error("covariate_names supplied as a single vector but components have different covariate counts (list them); supply a length-K list instead.")
    end
    isempty(supervised_lens) && return (names, intercept_flag)

    L = supervised_lens[1]

    if !intercept_specified && length(covariate_names) == L && length(covariate_names) != L - 1
        intercept_flag = false
        @warn "covariate_names lengths match gamma lengths exactly, implying no intercept; proceeding with intercept=false. Pass intercept=false explicitly to silence this warning."
    end

    expected = intercept_flag ? L - 1 : L
    length(covariate_names) == expected ||
        error("covariate_names must have length $expected; got $(length(covariate_names)).")

    cn = collect(String, covariate_names)
    for k in 1:K
        supervised[k] && (names[k] = cn)
    end

    return (names, intercept_flag)
end

function _dash_component_covariate_names(k, gamma_lens, intercept_flag, resolved_names)
    q = intercept_flag ? gamma_lens[k] - 1 : gamma_lens[k]
    q <= 0 && return String[]
    resolved_names[k] !== nothing && return resolved_names[k]
    return ["V$i" for i in 1:q]
end

function _dash_all_pairs_disjoint(name_sets)
    n = length(name_sets)
    for i in 1:n, j in (i+1):n
        isempty(intersect(name_sets[i], name_sets[j])) || return false
    end
    return true
end

function _dash_build_gamma_table(gamma_list, gamma_lens, supervised, resolved_names, K, intercept_flag)
    comp_names = [_dash_component_covariate_names(k, gamma_lens, intercept_flag, resolved_names) for k in 1:K]
    supervised_idx = findall(supervised)
    mode = _dash_all_pairs_disjoint([comp_names[k] for k in supervised_idx]) ? :pairs : :union

    rank = String[]
    note = String[]
    intercept_col = intercept_flag ? Union{Missing,Float64}[] : nothing

    if mode === :union
        gamma_col_names = String[]
        for k in supervised_idx, nm in comp_names[k]
            nm in gamma_col_names || push!(gamma_col_names, nm)
        end
        value_cols = Dict(nm => Union{Missing,Float64}[] for nm in gamma_col_names)

        for k in 1:K
            push!(rank, "R$k")
            if supervised[k]
                push!(note, "")
                covariate_vals = intercept_flag ? gamma_list[k][2:end] : gamma_list[k]
                intercept_flag && push!(intercept_col, gamma_list[k][1])
                row_vals = Dict{String,Any}(nm => missing for nm in gamma_col_names)
                for (i, nm) in enumerate(comp_names[k])
                    row_vals[nm] = covariate_vals[i]
                end
                for nm in gamma_col_names
                    push!(value_cols[nm], row_vals[nm])
                end
            else
                push!(note, "no covariates for this component")
                intercept_flag && push!(intercept_col, missing)
                for nm in gamma_col_names
                    push!(value_cols[nm], missing)
                end
            end
        end

        df = DataFrame(Rank=rank, Note=note)
        intercept_flag && (df.Intercept = intercept_col)
        for nm in gamma_col_names
            df[!, Symbol(nm)] = value_cols[nm]
        end
        return (df, :union)
    else
        covariates_col = String[]
        for k in 1:K
            push!(rank, "R$k")
            if supervised[k]
                push!(note, "")
                covariate_vals = intercept_flag ? gamma_list[k][2:end] : gamma_list[k]
                intercept_flag && push!(intercept_col, gamma_list[k][1])
                pieces = ["$(nm): $(round(covariate_vals[i]; digits=3))" for (i, nm) in enumerate(comp_names[k])]
                push!(covariates_col, join(pieces, ", "))
            else
                push!(note, "no covariates for this component")
                intercept_flag && push!(intercept_col, missing)
                push!(covariates_col, "")
            end
        end

        df = DataFrame(Rank=rank, Note=note)
        intercept_flag && (df.Intercept = intercept_col)
        df.Covariates = covariates_col
        return (df, :pairs)
    end
end

"""Resolve dash() inputs (fit or raw factors, labels, covariate names) into a static gamma table."""
function dash_data(; fit=nothing, Ef=nothing, Ew=nothing, lambda=nothing,
                    gamma_list=nothing, channel_names=nothing, time_labels=nothing,
                    channel_groups=nothing, covariate_names=nothing,
                    intercept::Union{Bool,Nothing}=nothing)
    intercept_specified = intercept !== nothing
    intercept_flag = intercept === nothing ? true : intercept

    raw_supplied = Ef !== nothing || Ew !== nothing || lambda !== nothing

    local K
    if raw_supplied
        (Ef !== nothing && Ew !== nothing && lambda !== nothing) ||
            error("When bypassing fit, all of Ef, Ew, and lambda must be supplied together.")
        Ef isa AbstractMatrix || error("Ef must be a matrix.")
        Ew isa AbstractMatrix || error("Ew must be a matrix.")
        K = size(Ef, 2)
        size(Ew, 2) == K ||
            error("Ew must have the same number of columns (K) as Ef; got $(size(Ew, 2)) vs $K.")
        length(lambda) == K ||
            error("lambda must have length K=$K; got $(length(lambda)).")
        if gamma_list !== nothing
            length(gamma_list) == K ||
                error("gamma_list must have length K=$K; got $(length(gamma_list)).")
        end
    elseif fit !== nothing
        nf = normalize_factors(fit)
        Ef = nf.Ef
        Ew = nf.Ew
        lambda = nf.lambda
        gamma_list = nf.gamma_list
        K = size(Ef, 2)
    else
        error("Either fit or all of Ef, Ew, and lambda must be provided.")
    end

    if channel_names !== nothing
        length(channel_names) == size(Ew, 1) ||
            error("channel_names must have length matching size(Ew, 1)=$(size(Ew, 1)); got $(length(channel_names)).")
    end
    if time_labels !== nothing
        length(time_labels) == size(Ef, 1) ||
            error("time_labels must have length matching size(Ef, 1)=$(size(Ef, 1)); got $(length(time_labels)).")
    end
    if channel_groups !== nothing
        length(channel_groups) == size(Ew, 1) ||
            error("channel_groups must have length matching size(Ew, 1)=$(size(Ew, 1)); got $(length(channel_groups)).")
    end

    gamma_lens = gamma_list === nothing ? zeros(Int, K) : [_dash_gamma_len(gamma_list[k]) for k in 1:K]
    supervised = gamma_lens .> 0

    resolved_names, intercept_flag = _dash_resolve_covariate_names(
        covariate_names, gamma_lens, supervised, K, intercept_flag, intercept_specified)

    any_supervised = any(supervised)

    gamma_table, gamma_table_mode = if any_supervised
        _dash_build_gamma_table(gamma_list, gamma_lens, supervised, resolved_names, K, intercept_flag)
    else
        (nothing, :none)
    end

    return (Ef=Ef, Ew=Ew, lambda=lambda, gamma_list=gamma_list,
            channel_names=channel_names, time_labels=time_labels,
            channel_groups=channel_groups, covariate_names=resolved_names,
            intercept=intercept_flag, K=K, any_supervised=any_supervised,
            gamma_table=gamma_table, gamma_table_mode=gamma_table_mode)
end

function _dash_prepare(payload, PlutoMod; dir=mktempdir())
    payload_path = joinpath(dir, "cosparta_dash_payload.jls")
    serialize(payload_path, payload)

    cell_codes = [
        """begin
	sidebar = Div([
		Div("Components"; class="cospar-card-header"),
		Div([sel_bond]; class="cospar-card-body"),
	]; class="cospar-sidebar")

	panel_card(title, body; body_class="cospar-card-body") = Div([
		Div(title; class="cospar-card-header"),
		Div([body]; class=body_class),
	]; class="cospar-card")

	panels = Div([
		panel_card("Temporal pattern", p_time),
		panel_card("Channel pattern", p_chan),
		panel_card("Component weight", p_weight),
		panel_card("Covariate coefficients", gamma_view; body_class="cospar-card-body gamma"),
	]; class="cospar-panels")

	root = Div([sidebar, panels]; class="cospar-root")

	@htl(\"\"\"
<style>
main { max-width: none !important; }
.cospar-dash { font-family: sans-serif; }
.cospar-dash h1 { font-size: 1.4rem; margin: 0 0 1rem 0; }
.cospar-root { display: flex; gap: 1rem; align-items: flex-start; }
.cospar-sidebar { flex: 0 0 220px; border: 1px solid #d0d0d0; border-radius: 8px; overflow: hidden; background: white; }
.cospar-sidebar plj-multi-checkbox { display: flex !important; flex-direction: column !important; }
.cospar-panels { flex: 1; display: grid; grid-template-columns: 1fr 1fr; grid-template-rows: 1fr 1fr; gap: 1rem; }
.cospar-card { border: 1px solid #d0d0d0; border-radius: 8px; overflow: hidden; background: white; }
.cospar-card-header { background: #f2f4f7; padding: 0.4rem 0.75rem; font-weight: 600; border-bottom: 1px solid #d0d0d0; }
.cospar-card-body { padding: 0.5rem; }
.cospar-card-body.gamma { overflow-x: auto; }
.cospar-card-body img { max-width: 100%; height: auto; }
</style>
<div class="cospar-dash">
	<h1>CoSparTA Fit Explorer</h1>
	\$root
</div>
\"\"\")
end""",
        "using PlutoUI, PlutoUI.ExperimentalLayout, Plots, DataFrames, Serialization, HypertextLiteral",
        """begin
	payload = deserialize($(repr(payload_path)))
	nothing
end""",
        """begin
	sel_bond = @bind sel_raw MultiCheckBox([k => "R\$k" for k in 1:payload.K]; default = collect(1:payload.K), select_all = true)
	nothing
end""",
        """begin
	sel = isempty(sel_raw) ? collect(1:payload.K) : sort(sel_raw)
	nothing
end""",
        """begin
	comp_palette = palette(:tab10)
	comp_color(k) = comp_palette[mod1(k, 10)]
	nothing
end""",
        """begin
	p_time = let
		tl = payload.time_labels === nothing ? collect(1:size(payload.Ef,1)) : payload.time_labels
		plt = plot(layout=(length(sel),1), legend=false, size=(640, 400))
		for (i,k) in enumerate(sel)
			plot!(plt, tl, payload.Ef[:,k]; sp=i, title="R\$k", seriestype=:path, color=comp_color(k),
				  xlabel=(i==length(sel) ? "Time" : ""), ylabel="Loading")
		end
		plt
	end
	nothing
end""",
        """begin
	p_chan = let
		cl = payload.channel_names === nothing ? string.(1:size(payload.Ew,1)) : string.(payload.channel_names)
		plt = plot(layout=(length(sel),1), legend=false, size=(640, 400))
		for (i,k) in enumerate(sel)
			bar!(plt, cl, payload.Ew[:,k]; sp=i, title="R\$k", color=comp_color(k),
				 xlabel=(i==length(sel) ? "Channels" : ""), ylabel="Loading", xrotation=45)
		end
		plt
	end
	nothing
end""",
        """begin
	p_weight = plot(["R\$k" for k in sel], payload.lambda[sel]; seriestype=:line, linecolor=:grey,
					 marker=:circle, markercolor=[comp_color(k) for k in sel], markerstrokecolor=[comp_color(k) for k in sel],
					 legend=false, xlabel="Component", ylabel="Weight", size=(640, 400))
	nothing
end""",
        """begin
	gamma_view = if !payload.any_supervised || payload.gamma_table === nothing
		md"This fit is fully unsupervised — no covariate coefficients to display."
	else
		wanted = Set("R\$k" for k in sel)
		payload.gamma_table[in.(payload.gamma_table.Rank, Ref(wanted)), :]
	end
	nothing
end""",
    ]

    cells = [Base.invokelatest(PlutoMod.Cell; code=code, code_folded=true) for code in cell_codes]
    nb = Base.invokelatest(PlutoMod.Notebook, cells)

    nb_path = joinpath(dir, "cosparta_dashboard.jl")
    nb.path = nb_path
    try
        Base.invokelatest(PlutoMod.save_notebook, nb)
    catch e
        e isa MethodError || rethrow()
        Base.invokelatest(PlutoMod.save_notebook, nb, nb_path)
    end

    return nb_path
end

function dash(fit::CoSparTAFit; kwargs...)
    return dash(; fit=fit, kwargs...)
end

"""Launch an interactive Pluto dashboard rendering the resolved dash_data() payload."""
function dash(; fit=nothing, Ef=nothing, Ew=nothing, lambda=nothing, gamma_list=nothing,
              channel_names=nothing, time_labels=nothing, channel_groups=nothing,
              covariate_names=nothing, intercept=nothing)
    payload = dash_data(fit=fit, Ef=Ef, Ew=Ew, lambda=lambda, gamma_list=gamma_list,
                         channel_names=channel_names, time_labels=time_labels,
                         channel_groups=channel_groups, covariate_names=covariate_names,
                         intercept=intercept)

    Base.find_package("Pluto") === nothing &&
        error("dash() requires the Pluto package. Install it with: import Pkg; Pkg.add(\"Pluto\"). The core package and dash_data() work without it.")

    PlutoMod = Base.require(Main, :Pluto)
    nb = _dash_prepare(payload, PlutoMod)

    println("Launching the CoSparTA dashboard in Pluto; this blocks the current session until the server is stopped.")
    Base.invokelatest(PlutoMod.run; notebook=nb)

    return nb
end

const _DASH_STIPPLE_APP_TEMPLATE = """
module CoSparTADashApp
using GenieFramework, PlotlyBase, StipplePlotly, Serialization
@genietools

const PAYLOAD = deserialize({{PAYLOAD_PATH}})
const NK = PAYLOAD.K
_ranks(sel) = isempty(sel) ? collect(1:NK) : sort(Int.(sel))

function weight_traces(sel)
    rs = _ranks(sel)
    [PlotlyBase.scatter(x = ["R\$k" for k in rs], y = PAYLOAD.lambda[rs],
                        mode = "lines+markers", name = "weight")]
end

@app begin
    @in selected = collect(1:NK)
    @out wtraces = weight_traces(collect(1:NK))
    @out wlayout = PlotlyBase.Layout(title = "Component weights")
    @onchange selected begin
        wtraces = weight_traces(selected)
    end
end

function ui()
    row([
        cell(class = "col-3", [
            card(class = "q-pa-md", [
                h6("Components")
                select(:selected,
                       options = [Dict("label" => "R\$k", "value" => k) for k in 1:NK],
                       multiple = true, emitvalue = true, mapoptions = true,
                       label = "Select components")
            ])
        ]),
        cell(class = "col-9", [
            card(class = "q-pa-md", [ plot(:wtraces, layout = :wlayout) ])
        ])
    ])
end

@page("/", ui)
Server.isrunning() || Server.up(async = true)
end
"""

"""Build the Genie/Stipple app-module file for the resolved dash_data() payload; does not launch it."""
function _dash_stipple_prepare(payload; dir=mktempdir())
    payload_path = joinpath(dir, "cosparta_dash_payload.jls")
    serialize(payload_path, payload)

    app_text = replace(_DASH_STIPPLE_APP_TEMPLATE,
                        "{{PAYLOAD_PATH}}" => repr(abspath(payload_path)))

    app_path = joinpath(dir, "cosparta_stipple_app.jl")
    write(app_path, app_text)

    return app_path
end

function dash_stipple(fit::CoSparTAFit; kwargs...)
    return dash_stipple(; fit=fit, kwargs...)
end

"""Launch an interactive Genie/Stipple dashboard rendering the resolved dash_data() payload."""
function dash_stipple(; fit=nothing, Ef=nothing, Ew=nothing, lambda=nothing, gamma_list=nothing,
                       channel_names=nothing, time_labels=nothing, channel_groups=nothing,
                       covariate_names=nothing, intercept=nothing)
    payload = dash_data(fit=fit, Ef=Ef, Ew=Ew, lambda=lambda, gamma_list=gamma_list,
                         channel_names=channel_names, time_labels=time_labels,
                         channel_groups=channel_groups, covariate_names=covariate_names,
                         intercept=intercept)

    Base.find_package("GenieFramework") === nothing &&
        error("dash_stipple() requires the GenieFramework, PlotlyBase, and StipplePlotly packages. Install them with: import Pkg; Pkg.add([\"GenieFramework\", \"PlotlyBase\", \"StipplePlotly\"]). The core package and dash_data() work without them.")

    app = _dash_stipple_prepare(payload)

    @info "Launching the CoSparTA Stipple dashboard at http://localhost:8000; to stop it, run `CoSparTADashApp.Server.down()`."
    Base.include(Main, app)

    return app
end
