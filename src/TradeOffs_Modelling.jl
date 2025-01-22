include("packages.jl")
include("TradeOffs_CommonCode.jl")

##Beverton Holt Model
function calc_m(para)
    @unpack a,b,K,p,τ = para
    m = (a-b*exp(-K*(τ+1)))*(p^(τ+1))
    return m
end

#Density independent survival of immature individuals
function BevertonHolt_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
        return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

let 
    time = 500
    timeseries = model_recursion(0.1,time,BevHoltPar(τ=3,p=0.55), BevertonHolt_model)
    plot(0:1:time,timeseries)
    # return timeseries
end

#issue with p=0.55 and τ=1 - shoots off to infinity

let     
    τrange = 2:1:7
    BHbifurcdata_p55 = model_τbifurc(τrange, BevertonHolt_model, BevHoltPar(), 0.55)
    BHbifurcdata_p5 = model_τbifurc(τrange, BevertonHolt_model, BevHoltPar(), 0.5)
    BHbifurcdata_p45 = model_τbifurc(τrange, BevertonHolt_model, BevHoltPar(), 0.45)
    scatter(τrange,BHbifurcdata_p55, label="p=0.55")
    scatter!(τrange,BHbifurcdata_p5, label="p=0.5")
    scatter!(τrange,BHbifurcdata_p45, label="p=0.45")
    xlabel!("τ")
    ylabel!("N")
end

#Code max 0 or equilibrium point for changing tau

#Cohort dependent survival of immature individuals (immature individuals exposed to density effects within cohort - but not density effect with mature individuals)
function BevertonHolt_modelII(Ndata, t, para)
    @unpack α,β,a,b,K,p,D,C,τ = para
        return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (D*((a-b*exp(-K*(τ+1)))*Ndata[t-τ]))/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*(a-b*exp(-K*(τ+1)))*Ndata[t-τ])
end

let 
    time = 500
    timeseries = model_recursion(0.1,time,BevHoltPar(τ=3,p=0.55), BevertonHolt_modelII)
    plot(0:1:time,timeseries)
    # return timeseries
end

#Mature dependent survival of immature individuals (immature individuals exposed to density effects with mature individuals)
#TODO

##Ricker Model
function Ricker_model(Ndata, t, para)
    @unpack r,β, a,b,K,p,τ = para
        return (Ndata[t] * exp(r*(1-(Ndata[t]/β)))) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

let 
    time = 500
    timeseries = model_recursion(0.1,time, RickerPar(τ=2,p=0.55), Ricker_model)
    plot(0:1:time,timeseries)
end

let     
    τrange = 0:1:10
    Rbifurcdata_p55 = model_τbifurc(τrange, Ricker_model, RickerPar(), 0.55)
    Rbifurcdata_p5 = model_τbifurc(τrange, Ricker_model, RickerPar(), 0.5)
    Rbifurcdata_p45 = model_τbifurc(τrange, Ricker_model, RickerPar(), 0.45)
    scatter(τrange, Rbifurcdata_p55, label="p=0.55")
    scatter!(τrange, Rbifurcdata_p5, label="p=0.5")
    scatter!(τrange, Rbifurcdata_p45, label="p=0.45")
    xlabel!("τ")
    ylabel!("N")
end

#Cohort dependent survival of immature individuals (immature individuals exposed to density effects within cohort - but not density effect with mature individuals)
#Attempt where mix Ricker model with bevertonholt solution to fraction of immature individuals that survive to t+1
function Ricker_modelIIa(Ndata, t, para)
    @unpack r,β,a,b,K,p,D,C,τ = para
        return (Ndata[t] * exp(r*(1-(Ndata[t]/β)))) + (D*((a-b*exp(-K*(τ+1)))*Ndata[t-τ]))/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*(a-b*exp(-K*(τ+1)))*Ndata[t-τ])
end

#TODO change ricker to same formulation as streipert paper -equation 30 without the births?

let 
    time = 500
    timeseries = model_recursion(0.1,time, RickerPar(τ=2,p=0.55), Ricker_modelIIa)
    plot(0:1:time,timeseries)
end

#Mature dependent survival of immature individuals (immature individuals exposed to density effects with mature individuals)
#TODO