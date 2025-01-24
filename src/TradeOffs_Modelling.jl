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
end

let
    τrange = 2:1:7
    BHIorbitdata = model_τorbit(τrange, BevertonHolt_model, BevHoltPar(p=0.55), 50)
    BHIIorbitdata = model_τorbit(τrange, BevertonHolt_modelII, BevHoltPar(), 50)
    test = plot(ylims=(0,5), xlims=(2,7))
    plot_combination(τrange, BHIorbitdata,:red)
    plot_combination(τrange, BHIIorbitdata,:black)
    xlabel!("τ")
    ylabel!("N")
end

#Mature dependent survival of immature individuals (immature individuals exposed to density effects with mature individuals)
#TODO

##Ricker Model
function Ricker_model(Ndata, t, para)
    @unpack d,c, a,b,K,p,τ = para
        return (Ndata[t] * exp(-d-c*Ndata[t])) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

let 
    time = 500
    timeseries = model_recursion(0.1,time, RickerPar(τ=3,p=0.55,c=0.05,d=0.05), Ricker_model)
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
    @unpack d,c,a,b,K,p,D,C,τ = para
        return (Ndata[t] * exp(-d-c*Ndata[t])) + (D*((a-b*exp(-K*(τ+1)))*Ndata[t-τ]))/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*(a-b*exp(-K*(τ+1)))*Ndata[t-τ])
end


let 
    time = 500
    timeseries = model_recursion(0.1,time, RickerPar(τ=2,p=0.55), Ricker_modelIIa)
    plot(0:1:time,timeseries)
end

let
    τrange = 0:1:10
    RIorbitdata = model_τorbit(τrange, Ricker_model, RickerPar(), 50)
    RIIorbitdata = model_τorbit(τrange, Ricker_modelIIa, RickerPar(), 50)
    test = plot(ylims=(0,5), xlims=(0,10))
    plot_combination(τrange, RIorbitdata,:red)
    plot_combination(τrange, RIIorbitdata,:black)
    xlabel!("τ")
    ylabel!("N")
end


#Are we expecting oscillations in Ricker or Beverton Holt model? and with what parameter values?

#Leslie matrix version of the Ricker model
function adultsurvival(A, para)
    @unpack d,c = para
    return exp(-d-c*A)
end

function juvenilebirth(A, para)
    @unpack D,C,a,b,K,p,τ = para
    R = a-b*exp(-K*(τ+1))
    return R*exp(-D-C*R*A)
end

function juvenilesurvival(J, para)
    @unpack D,C = para
    return exp(-D-C*J)
end


function tau3matrix(vector, para)
    return [[adultsurvival(vector[1], para), juvenilebirth(vector[1], para), 0, 0] [0,0,juvenilesurvival(vector[2], para),0]   [0,0,0,juvenilesurvival(vector[3], para)]   [juvenilesurvival(vector[4], para),0,0,0]]
end 

function model_Leslierecursion(vector0, time, para, Lesliematrix)
    AJvector = [Vector{Float64}() for _ in 1:time+1]
    AJvector[1] = vector0
    for t in 1:time
        AJvector[t+1] = Lesliematrix(AJvector[t], para)*AJvector[t]
    end
    return AJvector
end

function first_elements(vec_of_vecs)
    return [vec[1] for vec in vec_of_vecs]
end

let 
    time = 500
    timeseries = model_Leslierecursion([0.1,0.1,0.1,0.1], time, RickerPar(τ=3), tau3matrix)
    plot(0:1:time,first_elements(timeseries))
end

#Mature dependent survival of immature individuals (immature individuals exposed to density effects with mature individuals)
#TODO