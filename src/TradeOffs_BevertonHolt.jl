include("packages.jl")

@with_kw mutable struct BevHoltPar
    α::Float64 = 0.1
    β::Float64 = 0.3
    a::Float64 = 10
    b::Float64 = 200
    K::Float64 = 1
    p::Float64 = 0.4
    τ::Int = 5
end

function calc_m(para)
    @unpack a,b,K,p,τ = para
    m = (a-b*exp(-K*(τ+1)))*(p^(τ+1))
    return m
end

calc_m(BevHoltPar())

function BevertonHolt_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
        return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

function BevertonHolt_recursion(N0,time,para)
    @unpack τ = para
    N = fill(N0, τ+1)
    for t in τ+1:τ+time
        newN = BevertonHolt_model(N,t,para)
        append!(N,newN)
    end
    return N[τ+1:end]
end

let 
    time = 500
    timeseries = BevertonHolt_recursion(0.1,time,BevHoltPar(τ=2,p=0.55))
    plot(0:1:time,timeseries)
end