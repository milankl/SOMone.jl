# SOMone.jl
A 1D self-organizing map. Caution: This repository is in a preliminary state and presumably does not yield expected results. Feel free to improve it.

# Usage

Given a 1D data array, you specify the number of nodes `n` and the number of iterative steps `smax`, then
```
data = randn(500)   
n = 20
smax = 200
nodes = SOM(n,smax,data)
```
will return an array `nodes`, which are the locations of the `n` nodes.

# Problems

Current problems are: Agglomeration (several nodes tend to collapse to a single one), the SOM does not fit the data properly (without the maximum entropy initilisation the map is not able to shift its nodes onto a data set that is either offset, or has a much larger/smaller spread, the SOM doesn't converge (estimating the average rounding error on every step does not necessarily yield a smaller error for more SOM steps). 

# Why does it still work somewhat?

Most of the quantization of the data is done with a maximum entropy initilisation (src/initialisation.jl). At the moment the SOM just moves around these nodes a little bit.

