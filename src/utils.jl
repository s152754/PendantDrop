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