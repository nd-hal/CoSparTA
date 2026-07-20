struct TimeFactorPlot{T,L}
    values::T
    labels::L
    ranks::Vector{Int}
end

@recipe function f(plot::TimeFactorPlot)
    layout := (length(plot.ranks), 1)
    legend := false
    xlabel := "Time"
    ylabel := "Loading"
    for (panel, rank) in enumerate(plot.ranks)
        @series begin
            subplot := panel
            seriestype := :path
            title := "R$rank"
            plot.labels, plot.values[:, panel]
        end
    end
end

struct ChannelFactorPlot{T,L,G}
    values::T
    labels::L
    groups::G
    ranks::Vector{Int}
end

@recipe function f(plot::ChannelFactorPlot)
    layout := (length(plot.ranks), 1)
    legend := false
    xlabel := "Channels"
    ylabel := "Loading"
    for (panel, rank) in enumerate(plot.ranks)
        @series begin
            subplot := panel
            seriestype := :bar
            title := "R$rank"
            plot.labels, plot.values[:, panel]
        end
    end
end

"""Create a RecipesBase-compatible faceted view of time-mode factors."""
function plot_time_factors(result::Union{Nothing,CoSparTAFit}=nothing; ranks=nothing,
                           time_labels=nothing, normalize=true, Ef=nothing, kwargs...)
    Ef === nothing && result === nothing && throw(ArgumentError("provide a fit or Ef"))
    values = Ef === nothing ? (normalize ? result.U2_normed.mean : result.U2.mean) : Matrix(Ef)
    selected = ranks === nothing ? collect(axes(values, 2)) : Int.(ranks)
    labels = time_labels === nothing ? collect(axes(values, 1)) : collect(time_labels)
    length(labels) == size(values, 1) || throw(DimensionMismatch("time_labels must match Ef rows"))
    TimeFactorPlot(values[:, selected], labels, selected)
end

"""Create a RecipesBase-compatible faceted view of channel-mode factors."""
function plot_channel_factors(result::Union{Nothing,CoSparTAFit}=nothing; ranks=nothing,
                              channel_names=nothing, channel_groups=nothing,
                              normalize=true, Ew=nothing, kwargs...)
    Ew === nothing && result === nothing && throw(ArgumentError("provide a fit or Ew"))
    values = Ew === nothing ? (normalize ? result.U3_normed.mean : result.U3.mean) : Matrix(Ew)
    selected = ranks === nothing ? collect(axes(values, 2)) : Int.(ranks)
    labels = channel_names === nothing ? collect(axes(values, 1)) : collect(channel_names)
    length(labels) == size(values, 1) || throw(DimensionMismatch("channel_names must match Ew rows"))
    groups = channel_groups === nothing ? nothing : collect(channel_groups)
    groups !== nothing && length(groups) != size(values, 1) && throw(DimensionMismatch("channel_groups must match Ew rows"))
    ChannelFactorPlot(values[:, selected], labels, groups, selected)
end
