include("packages.jl")
include("TradeOffs_CommonCode.jl")

##Beverton Holt Model
#Density independent survival of immature individuals
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

# function timeembedding(tauval, Ndata)

let 
    time = 5000
    finalts=4900
    timeseries = model_recursion(0.1,time,RickerPar(τ=3,p=0.6, a=11.0), RickerConstant_model)
    plot(timeseries[finalts-3:end-3], timeseries[finalts:end])
    xlabel!("N(t-3)")
    ylabel!("N(t)")
    # return timeseries[4502:5001]
end


let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=2), 50))
    tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=3), 50))
    tau4data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=4), 50))
    tau5data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=5), 50))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,250.0)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,250.0)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,250.0)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,250.0)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
end

let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=2,p=0.6), 50))
    tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=3,p=0.6), 50))
    tau4data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=4,p=0.6), 50))
    tau5data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=5,p=0.6), 50))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,20.0)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,20.0)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,20.0)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,20.0)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
end

let
    τrange = 2:1:10
    RIorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", RickerPar(a=11.0,p=0.6), 150; upperval=10))
    scatter(RIorbitdata[1], RIorbitdata[2],color=:red, label="Ricker (constant): a=11.0,p=0.6")
    xlabel!("τ")
    ylabel!("N")
end

let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "p", RickerPar(τ=2, a=15.0), 150))
    tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "p", RickerPar(τ=3, a=15.0), 150))
    tau4data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "p", RickerPar(τ=4, a=15.0), 150))
    tau5data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "p", RickerPar(τ=5, a=15.0), 150))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("p")
    ylabel!("N")
    xlims!(0.0,1.0)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("p")
    ylabel!("N")
    xlims!(0.0,1.0)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("p")
    ylabel!("N")
    xlims!(0.0,1.0)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("p")
    ylabel!("N")
    xlims!(0.0,1.0)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
end

let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "α", RickerPar(τ=2, a=25.0), 150; upperval=2.0))
    println(calc_m(RickerPar(τ=2, a=25.0)))
    tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "α", RickerPar(τ=3, a=25.0), 150; upperval=2.0))
    println(calc_m(RickerPar(τ=3, a=25.0)))
    tau4data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "α", RickerPar(τ=4, a=25.0), 150; upperval=2.0))
    println(calc_m(RickerPar(τ=4, a=25.0)))
    tau5data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "α", RickerPar(τ=5, a=25.0), 150; upperval=2.0))
    println(calc_m(RickerPar(τ=5, a=25.0)))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("α")
    ylabel!("N")
    xlims!(0.0,2.0)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("α")
    ylabel!("N")
    xlims!(0.0,2.0)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("α")
    ylabel!("N")
    xlims!(0.0,2.0)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("α")
    ylabel!("N")
    xlims!(0.0,2.0)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
end

let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "β", RickerPar(τ=2, a=25.0), 50; upperval=0.5))
    println(calc_m(RickerPar(τ=2, a=25.0)))
    tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "β", RickerPar(τ=3, a=25.0), 50; upperval=0.5))
    println(calc_m(RickerPar(τ=3, a=25.0)))
    tau4data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "β", RickerPar(τ=4, a=25.0), 50; upperval=0.5))
    println(calc_m(RickerPar(τ=4, a=25.0)))
    tau5data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "β", RickerPar(τ=5, a=25.0), 50; upperval=0.5))
    println(calc_m(RickerPar(τ=5, a=25.0)))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("β")
    ylabel!("N")
    xlims!(0.0,0.5)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("β")
    ylabel!("N")
    xlims!(0.0,0.5)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("β")
    ylabel!("N")
    xlims!(0.0,0.5)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("β")
    ylabel!("N")
    xlims!(0.0,0.5)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
end


let     
    τrange = 0:1:10
    Rbifurcdata_p55 = model_τbifurc(τrange, RickerConstant_model, RickerPar(), 0.55)
    Rbifurcdata_p5 = model_τbifurc(τrange, RickerConstant_model, RickerPar(), 0.5)
    Rbifurcdata_p45 = model_τbifurc(τrange, RickerConstant_model, RickerPar(), 0.45)
    scatter(τrange, Rbifurcdata_p55, label="p=0.55")
    scatter!(τrange, Rbifurcdata_p5, label="p=0.5")
    scatter!(τrange, Rbifurcdata_p45, label="p=0.45")
    xlabel!("τ")
    ylabel!("N")
end

#Cohort dependent survival of immature individuals (immature individuals exposed to density effects within cohort - but not density effect with mature individuals)
#Attempt where mix Ricker model with bevertonholt solution to fraction of immature individuals that survive to t+1

let 
    time = 500
    timeseries = model_recursion(0.1,time, RickerPar(τ=2.0, a=2.0, b=1.0, p=0.8, d=0.0001, K=1.0, c=0.3), RickerBeverton_model)
    plot(0:1:time,timeseries)
end

let
    τrange = 0:1:10
    RIorbitdata = model_τorbit(τrange, RickerConstant_model, RickerPar(), 50)
    RIIorbitdata = model_τorbit(τrange, RickerBeverton_model, RickerPar(), 50)
    test = plot(ylims=(0,5), xlims=(0,10))
    plot_combination(τrange, RIorbitdata,:red)
    plot_combination(τrange, RIIorbitdata,:black)
    xlabel!("τ")
    ylabel!("N")
end

#QUESTION - what is the constraint for RickerBeverton model - already constrained by beverton part of density dependence. EXCEPT NEED LOWER BOUND TO a so that g(τ) does not go negative (negative g makes no biological sense)
let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=2), 50; upperval=50))
    tau3data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=3), 50; upperval=50))
    tau4data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=4), 50; upperval=50))
    tau5data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=5), 50; upperval=50))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,50.0)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,50.0)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,50.0)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,50.0)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
end

#Leslie matrix version of the Ricker model
function adultsurvival(A, para)
    @unpack α,β = para
    return exp(-α-β*A)
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
    return [τrange,data]
end

let 
    timeseries = model_recursion(0.1,5000,RickerPar(τ=0, a=300.0,α=1.5,β=2.0,D=0.5,C=0.15), RickerLeslie_τ0_model)
    return timeseries[end-100:end]
end

let
    τrange = 1:1:15
    # RIorbitdata = model_τorbit(τrange, RickerConstant_model, RickerPar(a=100, p=0.6), 50)
    # RIIorbitdata = model_τorbit(τrange, RickerBeverton_model, RickerPar(a=100,D=0.5,C=0.15), 50)
    RLorbitdata = flattenorbitdata(LeslieMatrixOrbitDiagram(τrange, 50000, 100, RickerPar(a=300.0,α=1.5,β=2.0,D=0.5,C=0.15), 0.1, LeslieMatrix))
    RLτ0data = model_recursion(0.1,5000,RickerPar(τ=0, a=300.0,α=1.5,β=2.0,D=0.5,C=0.15), RickerLeslie_τ0_model)
    RLτ0data_trans=RLτ0data[end-100:end]
    p2=scatter(RLorbitdata[1], RLorbitdata[2],color=:black)
    scatter!(zeros(length(RLτ0data_trans)), RLτ0data_trans, color=:black)
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0,15.0)
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
    RIorbitdata = model_τorbit(τrange, RickerConstant_model, RickerPar(p=0.6), 50)
    RIIorbitdata = model_τorbit(τrange, RickerBeverton_model, RickerPar(), 50)
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