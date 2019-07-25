""" Find the Best Matching Unit (BMU) given an array dist of distances."""
function bestMatchingUnit(dist::Array{Float64,1})
    return argmin(abs.(dist))
end

function neighbourhood(bmu::Int,s::Int,smax::Int,nodes::Array{Float64,1},w::Float64)

    decay = αexp(s,1e-3)
    shrink = αexp(s,5e-3)
    # decay = αlin(s,smax)
    # shrink = αlin(s,smax)

    # Gaussian neighbourhood function
    d = decay*exp.(-(((nodes .- nodes[bmu])/(w*shrink)).^2)./2)
    return d
end
