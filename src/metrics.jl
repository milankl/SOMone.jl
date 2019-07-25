# """ Sign preserving 1D distance between nodes and a data sample x. """
# function distance!(dist::Array{Float64,1},nodes::Array{Float64,1},x::Float64)
#     n = length(nodes)
#     @boundscheck n == length(dist) || throw(BoundsError())
#     @inbounds for i in 1:n
#         dist[i] = x - nodes[i]
#     end
# end

""" Sign preserving 1D distance between nodes and a data sample x. """
function distance(nodes::Array{Float64,1},x::Float64)
    return x .- nodes
end

""" Find the Best Matching Unit (BMU) given an array dist of distances."""
function bestMatchingUnit(dist::Array{Float64,1})
    return argmin(abs.(dist))
end

function αexp(s::Int,r::Float64=0.05)
    return exp(-r*s)
end

function αlin(s::Int,smax::Int)
    return 1-(s-1)/smax
end
