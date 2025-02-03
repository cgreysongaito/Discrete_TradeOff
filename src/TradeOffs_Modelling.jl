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
function BevHoltI_equi(para)
    @unpack α,β = para
    m = calc_m(para)
    equi = ((1+α)*m-α)/(β*(1-m))
    if equi>0
        return equi
    else
        return 0
    end
end

let 
    datap040=[BevHoltI_equi(BevHoltPar(τ=τval,p=0.4)) for τval in 0:1:6]
    datap050=[BevHoltI_equi(BevHoltPar(τ=τval,p=0.5)) for τval in 0:1:6]
    datap055=[BevHoltI_equi(BevHoltPar(τ=τval,p=0.55)) for τval in 0:1:6]
    equifig = scatter(0:1:6,datap040, label="p=0.4")
    scatter!(0:1:6,datap050, label="p=0.5")
    scatter!(0:1:6,datap055, label="p=0.55")
    xlabel!("τ")
    ylabel!("N*(τ)")
end


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
    g=a-(b*exp(-K*(τ+1)))
        return (Ndata[t] * exp(-d-c*Ndata[t])) + g*(p^(τ+1))*Ndata[t-τ]
end

let 
    time = 1000
    timeseries = model_recursion(50.0,time, RickerPar(τ=6,a=15,p=0.6,C=0.15,D=0.5), Ricker_model)
    plot(0:1:1000,timeseries)
    # plot(450:1:500,timeseries[450:501], color=:red)
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
    g=a-(b*exp(-K*(τ+1)))
        return (Ndata[t] * exp(-d-c*Ndata[t])) + (D/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*g*Ndata[t-τ]))*g*Ndata[t-τ]
end


let 
    time = 500
    timeseries = model_recursion(0.1,time, RickerPar(τ=2.0, a=2.0, b=1.0, p=0.8, d=0.0001, K=1.0, c=0.3), Ricker_modelIIa)
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

function adultdensity(A,para)
    @unpack D, C = para
    return exp(-D-C*A)
end

function LeslieMatrix(τval, para, AJvector)
    local_par = deepcopy(para)
    local_par.τ = τval
    matrix = zeros(Float64,τval+1,τval+1)
    matrix[1,1] = adultsurvival(AJvector[1], local_par)
    matrix[2,1] = juvenilebirth(AJvector[1], local_par)
    matrix[1,τval+1] = juvenilesurvival(AJvector[τval], local_par)
    for i in 3:τval+1
        matrix[i,i-1] = juvenilesurvival(AJvector[i-1], local_par)
    end
    return matrix
end

function model_Leslierecursion(τval, time, para, init, lesliematrix)
    initvector = fill(init, τval+1)
    AJvector = [Vector{Float64}() for _ in 1:time+1]
    AJvector[1] = initvector
    for t in 1:time
        AJvector[t+1] = lesliematrix(τval, para, AJvector[t])*AJvector[t]
    end
    return AJvector
end

function first_elements(vec_of_vecs)
    return [vec[1] for vec in vec_of_vecs]
end

let 
    time = 500
    τval=3
    timeseries = model_Leslierecursion(τval, time, RickerPar(τ=τval, a=400,D=0.5,C=0.15), 0.1, LeslieMatrix)
    plot(0:1:time,first_elements(timeseries))
    # return first_elements(timeseries)[end-50:end]
end


function LeslieMatrixOrbitDiagram(τrange, time, finalts, para, init, lesliematrix)
    data = Vector{Vector{Float64}}(undef, length(τrange))
    @threads for τi in eachindex(τrange)
        timeseries = model_Leslierecursion(τrange[τi], time, para, init, lesliematrix)
        data[τi] = first_elements(timeseries)[end-finalts:end]
    end
    return data
end

let
    τrange = 1:1:15
    # RIorbitdata = model_τorbit(τrange, Ricker_model, RickerPar(a=100, p=0.6), 50)
    # RIIorbitdata = model_τorbit(τrange, Ricker_modelIIa, RickerPar(a=100,D=0.5,C=0.15), 50)
    RLorbitdata = LeslieMatrixOrbitDiagram(τrange, 50000, 100, RickerPar(a=300,d=1.5,c=2.0,D=0.5,C=0.15), 0.1, LeslieMatrix)
    test = plot(ylims=(0.0,0.2), xlims=(0,15))
    # plot_combination(τrange, RIorbitdata,:black)
    # plot_combination(τrange, RIIorbitdata,:blue)
    plot_combination(τrange, RLorbitdata,:red)
    xlabel!("τ")
    ylabel!("N")
end

#Mature dependent survival of immature individuals (immature individuals exposed to density effects with mature individuals and their own cohort)
function LeslieMatrix_AdultCohort(τval, para, AJvector)
    local_par = deepcopy(para)
    local_par.τ = τval
    matrix = zeros(Float64,τval+1,τval+1)
    matrix[1,1] = adultsurvival(AJvector[1], local_par)
    matrix[2,1] = juvenilebirth(AJvector[1], local_par)
    matrix[1,τval+1] = adultdensity(AJvector[1], local_par)*juvenilesurvival(AJvector[τval], local_par)
    for i in 3:τval+1
        matrix[i,i-1] = adultdensity(AJvector[1], local_par)*juvenilesurvival(AJvector[i-1], local_par)
    end
    return matrix
end

let
    τrange = 2:1:15
    RIorbitdata = model_τorbit(τrange, Ricker_model, RickerPar(p=0.6), 50)
    RIIorbitdata = model_τorbit(τrange, Ricker_modelIIa, RickerPar(), 50)
    RLIorbitdata = LeslieMatrixOrbitDiagram(τrange, 500, 50, RickerPar(), 0.1, LeslieMatrix)
    RLIIorbitdata = LeslieMatrixOrbitDiagram(τrange, 500, 50, RickerPar(), 0.1, LeslieMatrix_AdultCohort)
    test = plot(ylims=(-0.1,8), xlims=(0,15))
    plot_combination(τrange, RIorbitdata,:black)
    plot_combination(τrange, RIIorbitdata,:blue)
    plot_combination(τrange, RLIorbitdata,:red)
    plot_combination(τrange, RLIIorbitdata,:orange)
    xlabel!("τ")
    ylabel!("N")
end