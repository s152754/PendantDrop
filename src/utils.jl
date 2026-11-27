using OrderedCollections
using Distributions

function dictsample(params::OrderedCollections.OrderedDict{String, Distribution}, nsample::Integer)
    Random.seed!(15)
    sample  = zeros(nsample, length(params))
    for (i, value) in enumerate(values(params))
        sample[:,i] = rand(value, nsample)
    end

    return sample
end

function dictsamplemean(means::Vector, stds::Vector, ncols::Int)
    Random.seed!(15)
    @assert length(means) == length(stds) "Means and stds must have the same length"
    nrows = length(means)
    mat = zeros(nrows, ncols)

    for i in 1:nrows
        dist = Normal(means[i], stds[i])
        mat[i, :] .= rand(dist, ncols)
    end

    return mat
end

function zoomplot_fn(zoom_factor, img)
    h, w = size(img)[1:2]   # height, width

    
    x_center = w / 2
    y_center = h / 2

    x_half = (w * zoom_factor) / 2
    y_half = (h * zoom_factor) / 2

    xlims = (x_center - x_half, x_center + x_half)
    ylims = (y_center - y_half, y_center + y_half)


    return xlims, ylims
end

function determine_arclength(rdfull::Vector{Int64}, zdfull::Vector{Int64})
    """
    Determine the arc length of a curve defined by points (rdfull, zdfull)

    Args:
        rdfull: Vector of radial coordinates
        zdfull: Vector of vertical coordinates
    Returns:
        ds: Vector of arc length increments between consecutive points
        sdfull: Vector of cumulative arc length values
    """
    n = length(rdfull)
    ds = [sqrt((rdfull[i+1] - rdfull[i])^2 + (zdfull[i+1] - zdfull[i])^2) for i in 1:n-1]
    sdfull = vcat(0.0, cumsum(ds))  # cumulative arc length starting at 0
    return ds, sdfull
end