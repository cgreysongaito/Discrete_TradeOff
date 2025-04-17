include("packages.jl")
include("TradeOffs_CommonCode.jl")

#TODO LIST
#make function that calculates lower and upperbounds of a for all tau
#make function that calculates lower and upper bounds of p for all tau

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
    using Random
    equifig = scatter(collect(0:1:6) .+ randn(7) .*0.05, datap040, label="p=0.4",markersize=6)
    scatter!(collect(0:1:6) .+ randn(7) .*0.05, datap050, label="p=0.5",markersize=6)
    scatter!(collect(0:1:6) .+ randn(7) .*0.05, datap055, label="p=0.55",markersize=6)
    xlabel!("τ")
    ylabel!("N*(τ)")
    plot!(grid=false)
    savefig(joinpath(abpath(), "figs/BevHoltI_equi.png"))
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

#BevertonHolt Ricker model
function BHadultsurvival(A, para)
    @unpack α,β = para
    return 1/(1+α+β*A)
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

function BRLeslieMatrix(τval, para, AJvector)
    local_par = deepcopy(para)
    local_par.τ = τval
    matrix = zeros(Float64,τval+1,τval+1)
    matrix[1,1] = BHadultsurvival(AJvector[1], local_par)
    matrix[2,1] = juvenilebirth(AJvector[1], local_par)
    matrix[1,τval+1] = juvenilesurvival(AJvector[τval], local_par)
    for i in 3:τval+1
        matrix[i,i-1] = juvenilesurvival(AJvector[i-1], local_par)
    end
    return matrix
end


BRLeslieMatrix(3, BevHoltPar(τ=3, a=25.0,α=1.5,β=2.0,D=0.5,C=0.15), [0.1, 0.2, 0.3, 0.4])
function BRmodel_Leslierecursion(τval, time, para, init, lesliematrix)
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
    time = 5000
    τval=1
    timeseries = BRmodel_Leslierecursion(τval, time, BevHoltPar(τ=τval, a=25.0,α=1.5,β=2.0,D=0.5,C=0.15), 0.1, BRLeslieMatrix)
    plot(0:1:time,first_elements(timeseries))
    # return first_elements(timeseries)[end-50:end]
end


function BRLeslieMatrixOrbitDiagram(τrange, time, finalts, para, init, lesliematrix)
    data = Vector{Vector{Float64}}(undef, length(τrange))
    @threads for τi in eachindex(τrange)
        timeseries = BRmodel_Leslierecursion(τrange[τi], time, para, init, lesliematrix)
        data[τi] = first_elements(timeseries)[end-finalts:end]
    end
    return [τrange,data]
end

function BevRickerLeslie_τ0_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(τ+1))
        return (Ndata[t] * 1/(1+α+β*Ndata[t])) + g*exp(-α-β*g*Ndata[t])*Ndata[t]
end

let
    τrange = 1:1:15
    # RIorbitdata = model_τorbit(τrange, RickerConstant_model, RickerPar(a=100, p=0.6), 50)
    # RIIorbitdata = model_τorbit(τrange, RickerBeverton_model, RickerPar(a=100,D=0.5,C=0.15), 50)
    par=BevHoltPar(τ=0, a=100.0,α=1.5,β=0.1,D=0.9,C=0.15)
    BRLorbitdata = flattenorbitdata(BRLeslieMatrixOrbitDiagram(τrange, 50000, 100, par, 0.1, BRLeslieMatrix))
    BRLτ0data = model_recursion(0.1,5000,par, BevRickerLeslie_τ0_model)
    BRLτ0data_trans=BRLτ0data[end-100:end]
    p2=scatter(BRLorbitdata[1], BRLorbitdata[2],color=:black, label="")
    scatter!(zeros(length(BRLτ0data_trans)), BRLτ0data_trans, color=:black, label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0,15.0)
    savefig(joinpath(abpath(), "figs/tauorbitdiagram_BHRicker.pdf"))
end

#Mature dependent survival of immature individuals (immature individuals exposed to density effects with mature individuals)

