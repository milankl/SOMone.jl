""" Sign preserving 1D distance between nodes and a data sample x. """
function distance(nodes::Array{Float64,1},x::Float64)
    return x .- nodes
end

function αexp(s::Int,r::Float64=0.05)
    return exp(-r*s)
end

function αlin(s::Int,smax::Int)
    return 1-(s-1)/smax
end

function maxspread(data::AbstractArray)
    return maximum(data)-minimum(data)
end

function crenel(x::Real, Δ::Real)
    if x > -Δ && x < Δ
        return one(x)
    else
        return zero(x)
    end
end


function average_rounding_error(data::Array{Float64,1},nodes::Array{Float64,1})
    n = length(data)
    S = 0.0
    for i = 1:n
        S += minimum(abs.(distance(nodes,data[i])))/n
    end
    return S
end
