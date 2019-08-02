module SOMone

export SOM,maxentropy

using Statistics, Printf, Random, StatsBase

include("metrics.jl")
include("neighbourhood.jl")
include("initialisation.jl")
include("selforganise.jl")

end
