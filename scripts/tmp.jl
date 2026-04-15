using ProgressMeter
using Base.Threads
using DrWatson
using Distributed

OUTDIR = "efficiency"

includet("tmp-pmax.jl")
@showprogress for id in 1:3
    runid(id)
end

@showprogress @threads for id in 1:200
    runid(id)
end