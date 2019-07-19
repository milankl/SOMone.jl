using PyPlot, Statistics, Printf


data = vcat(randn(500),randn(500).+4.0)

nnodes = 50

extr = extrema(data)
nodes = collect(range(extr[1],stop=extr[2],length=nnodes))

w = std(data)*5*length(nodes)/length(data)

""" Sign preserving 1D distance between nodes and a data sample x. """
function distance(nodes::Array{Float64,1},x::Float64)
    # n = length(nodes)
    # @boundscheck n == length(dist) || throw(BoundsError())
    # for i in 1:n
    #     @inbounds dist[i] = x - nodes[i]
    # end
    return x .- nodes
end

""" Find the Best Matching Unit (BMU) given an array dist of distances."""
function bestMatchingUnit(dist::Array{Float64,1})
    return argmin(abs.(dist))
end

function neighbourhood(bmu::Int,s::Int,nodes::Array{Float64,1},w::Float64)

    decay = α(s,5e-3)
    shrink = α(s,1e-1)

    # Gaussian neighbourhood function
    # if bmu == 1s
    #     w = abs(nodes[2]-nodes[1])
    # elseif bmu == length(nodes)
    #     w = abs(nodes[end]-nodes[end-1])
    # else
    #     w = abs(nodes[bmu+1]-nodes[bmu-1])/2
    # end

    d = decay*exp.(-(((nodes .- nodes[bmu])/(w*shrink)).^2)./2)
    return d
end

function α(s::Int,r::Float64=0.05)
    return exp(-r*s)
end


function update!(s::Int,data::Array{Float64,1},nodes::Array{Float64,1},w::Float64=1.0)
    x = rand(data)
    dist = distance(nodes,x)
    bmu = bestMatchingUnit(dist)
    θ = neighbourhood(bmu,s,nodes,w)
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

##

fig,ax = subplots(1,1,figsize=(8,2))
ax.plot(data,ones(length(data)),"|",lw=0.1)
l1, = ax.plot(nodes,zeros(length(nodes)),".-")
yticks([])
ylim(-1,2)
title("s = 0", loc="left")
tight_layout()
for s in 1:300
    update!(s,data,nodes,w)
    er = average_rounding_error(data,nodes)
    ser = @sprintf("%.4f",er)
    title("s = $s, r = $ser", loc="left")
    l1.set_data(nodes,zeros(length(nodes)))
    pause(α(s,0.5)+0.05)
end
