""" Find the Best Matching Unit (BMU) given an array dist of distances."""
function bestMatchingUnit(dist::Array{Float64,1})
    return argmin(abs.(dist))
end

function neighbourhood(bmu::Int,s::Int,smax::Int,nodes::Array{Float64,1},w::Float64,t::Float64,r::Float64)

    decay = αexp(s,t)
    shrink = αexp(s,r)

    # Gaussian neighbourhood function
    d = decay*exp.(-(((nodes .- nodes[bmu])/(w*shrink)).^2)./2)
    return d
end
