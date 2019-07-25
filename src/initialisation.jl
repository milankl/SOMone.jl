""" Preconditioner based on maximum entropy method. """
function maxentropy(r::Int,data::Array{Float64,1})

    N = length(data)
    n = N ÷ r

    # throw away random data for equally sized chunks of data
    data = shuffle(data[:])[1:n*r]
    sort!(data)

    nodes = Array{Float64,1}(undef,r)
    nodes[1] = minimum(data)
    nodes[2] = (2*data[1] + data[n] + data[n+1])/4

    for i in 2:r-1
        nodes[i+1] = (data[(i-1)*n] + data[(i-1)*n + 1] + data[i*n] + data[i*n + 1])/4
    end

    nodes[end-1] = (2*data[r*n] + data[(r-1)*n] + data[(r-1)*n-1])/4
    nodes[end] = maximum(data)

    return nodes
end
