function update!(s::Int,smax::Int,data::Array{Float64,1},nodes::Array{Float64,1},w::Float64,t::Float64,r::Float64)
    x = rand(data)
    dist = distance(nodes,x)
    bmu = bestMatchingUnit(dist)
    θ = neighbourhood(bmu,s,smax,nodes,w,t,r)
    for i in 1:length(nodes)
        nodes[i] += θ[i]*dist[i]
    end
end

function SOM(n::Int,smax::Int,data::Array{Float64,1};
            w::Float64=0.1,     # width scale of the Gaussian neighbourhood
            t::Float64=0.001,    # decay rate
            r::Float64=0.001)    # shrink rate

    nodes = maxentropy(n,data)

    w = maxspread(data)/length(nodes)*w

    for s in 1:smax
        update!(s,smax,data,nodes,w,t,r)
    end

    return nodes
end
