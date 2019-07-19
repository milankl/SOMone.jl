using PyPlot, Statistics, Printf, Random

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

function neighbourhood(bmu::Int,s::Int,smax::Int,nodes::Array{Float64,1},w::Float64)

    decay = αexp(s,1e-5)
    shrink = αexp(s,5e-5)
    # decay = αlin(s,smax)
    # shrink = αlin(s,smax)

    # Gaussian neighbourhood function
    d = decay*exp.(-(((nodes .- nodes[bmu])/(w*shrink)).^2)./2)
    return d
end

function αexp(s::Int,r::Float64=0.05)
    return exp(-r*s)
end

function αlin(s::Int,smax::Int)
    return 1-(s-1)/smax
end

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

data = vcat(randn(50000),randn(50000).+4.0)

Nnodes = 32767
nodes = maxentropy(Nnodes,data)

w = std(data)*length(nodes)/length(data)/40

smax = 10000
er = average_rounding_error(data,nodes)
ser = @sprintf("%.8f",er)
println("\nStart: $ser")
# global ermin = er

for s in 1:smax
    update!(s,smax,data,nodes,w)
    # er = average_rounding_error(data,nodes)
    #
    # global ermin
    # if er < ermin
    #     ermin = er
    # end
    #
    # ser = @sprintf("%.6f",er)
    print("\r\u1b[K")
    progress = Int(round(s/smax*100))
    print("$progress%")
end
er = average_rounding_error(data,nodes)
ser = @sprintf("%.8f",er)
println(", End: $ser")

##

fig,ax = subplots(1,1,figsize=(8,2))
ax.plot(data,ones(length(data)),"|",lw=0.1)
l1, = ax.plot(nodes,zeros(length(nodes)),".-")
yticks([])
ylim(-1,2)
title("s = 0", loc="left")
tight_layout()

# smax = 300
# for s in 1:smax
#     update!(s,smax,data,nodes,w)
#     er = average_rounding_error(data,nodes)
#     ser = @sprintf("%.4f",er)
#     title("s = $s, r = $ser", loc="left")
#     l1.set_data(nodes,zeros(length(nodes)))
#     pause(αexp(s,0.5)+0.05)
# end
