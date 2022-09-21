using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using Combinatorics, StatsBase

function bitarr_to_int5(arr,s=0)
    v = 1
    for i in view(arr,length(arr):-1:1)
        s += v*i
        v <<= 1
    end 
    s
end

function trans(x, v::Dict{T, Int}, l) where T
    z = collect(1:l)
    idxs = Vector{Int}[]
    for k in x
        push!(idxs, z[k])
        deleteat!(z, k)
    end
    res = Vector{T}(undef, l)
    for (j, k) in enumerate(keys(v))
        for i in idxs[j]
            res[i] = k
        end
    end
    res
end

function myperms(x)
    v = countmap(x)
    s = Int[length(x)]
    for (k,y) in v
        l = s[end]-y
        l > 0 && push!(s, l)
    end
    iter = Iterators.product((combinations(1:s[i], vv) for (i, vv) in enumerate(values(v)))...)
    (trans(z, v, length(x)) for z in iter)
end

num_ones = 13
num_zeros = 9

## CODE ## 
a = zeros(num_ones+num_zeros)
a[1:num_ones] .= 1 #[111000]
all_perms = myperms(a) # all unique permutations of a 
labels = bitarr_to_int5.(all_perms) # vectorized bit arrary -to- int 


@show labels 




