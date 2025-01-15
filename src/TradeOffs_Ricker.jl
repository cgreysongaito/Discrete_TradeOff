include("packages.jl")
@with_kw mutable struct RickerPar
    r::Float64 = 0.1
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

function Ricker_model(Ndata, t, para)
    @unpack r,β, a,b,K,p,τ = para
        return (Ndata[t] * exp(r*(1-(Ndata[t]/β)))) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

#Question is the k in the ricker model different from the K in the second half of this equation? I am assuming different

function Ricker_recursion(N0,time,para)
    @unpack τ = para
    N = fill(N0, τ+1)
    for t in τ+1:τ+time
        newN = Ricker_model(N,t,para)
        append!(N,newN)
    end
    return N[τ+1:end]
end

let 
    time = 500
    timeseries = Ricker_recursion(0.1,time,RickerPar(τ=2,p=0.55))
    plot(0:1:time,timeseries)
end

function Ricker_τbifurc(τrange, pval)
    data= zeros(length(τrange))
    @threads for i in eachindex(τrange)
        timeseries = Ricker_recursion(0.1,500,RickerPar(τ=τrange[i],p=pval))
        data[i] = timeseries[end]
    end
    return data
end

let     
    τrange = 0:1:10
    bifurcdata_p55 = Ricker_τbifurc(τrange,0.55)
    bifurcdata_p5 = Ricker_τbifurc(τrange,0.5)
    bifurcdata_p45 = Ricker_τbifurc(τrange,0.45)
    plot(τrange,bifurcdata_p55, label="p=0.55")
    plot!(τrange,bifurcdata_p5, label="p=0.5")
    plot!(τrange,bifurcdata_p45, label="p=0.45")
    xlabel!("τ")
    ylabel!("N")
end