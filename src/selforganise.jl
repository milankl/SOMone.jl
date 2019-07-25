function update!(s::Int,smax::Int,data::Array{Float64,1},nodes::Array{Float64,1},w::Float64=1.0)
    x = rand(data)
    dist = distance(nodes,x)
    bmu = bestMatchingUnit(dist)
    θ = neighbourhood(bmu,s,smax,nodes,w)
    for i in 1:length(nodes)
        nodes[i] += θ[i]*dist[i]
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
