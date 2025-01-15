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

#Question what is the history before t=0?

let 
    time = 500
    timeseries = BevertonHolt_recursion(0.1,time,BevHoltPar(τ=1,p=0.55))
    # plot(0:1:time,timeseries)
    return timeseries
end

function BevertonHolt_τbifurc(τrange, pval)
    data= zeros(length(τrange))
    @threads for i in eachindex(τrange)
        timeseries = BevertonHolt_recursion(0.1,500,BevHoltPar(τ=τrange[i],p=pval))
        data[i] = timeseries[end-50]
    end
    return data
end

let     
    τrange = 2:1:7
    bifurcdata_p55 = BevertonHolt_τbifurc(τrange,0.55)
    bifurcdata_p5 = BevertonHolt_τbifurc(τrange,0.5)
    bifurcdata_p45 = BevertonHolt_τbifurc(τrange,0.45)
    scatter(τrange,bifurcdata_p55, label="p=0.55")
    scatter!(τrange,bifurcdata_p5, label="p=0.5")
    scatter!(τrange,bifurcdata_p45, label="p=0.45")
    xlabel!("τ")
    ylabel!("N")
end

#issue with p=0.55 and τ=1 - shoots off to infinity