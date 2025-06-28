include("packages.jl")
include("TradeOffs_CommonCode.jl")
using LaTeXStrings
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
    datap040=[BevHoltI_equi(BevHoltPar(τ=τval,p=0.4)) for τval in 0:1:7]
    datap050=[BevHoltI_equi(BevHoltPar(τ=τval,p=0.5)) for τval in 0:1:7]
    datap055=[BevHoltI_equi(BevHoltPar(τ=τval,p=0.55)) for τval in 0:1:7]
    using Random
    equifig = scatter(collect(0:1:7) .+ randn(8) .*0.05, datap040, label="\${\\tilde{p}}=0.4\$",markersize=6, color="#440154FF")
    plot!(collect(0:1:7), datap040, label="",linestyle=:dash, color="#440154FF")
    scatter!(collect(0:1:7) .+ randn(8) .*0.05, datap050, label="\${\\tilde{p}}=0.5\$",markersize=6, color="#238A8DFF")
    plot!(collect(0:1:7), datap050, label="",linestyle=:dash, color="#238A8DFF")
    scatter!(collect(0:1:7) .+ randn(8) .*0.05, datap055, label="\${\\tilde{p}}=0.55\$",markersize=6, color="#FDE725FF")
    plot!(collect(0:1:7), datap055, label="",linestyle=:dash, color="#FDE725FF")
    xlabel!("τ")
    ylabel!("N*(τ)")
    savefig(joinpath(abpath(), "figs/BevHoltI_equi.pdf"))
end


#Cohort dependent survival of immature individuals (immature individuals exposed to density effects within cohort - but not density effect with mature individuals)

#root solve BevertonholtBevertonholt equality
function BHBH_existence_check(para)
    @unpack α,β,D,C,a,b,K,τ = para
    g = calc_g(para)
    return (1-(1/1+α))-(g/(1+D)^(τ+1))
end

function Nequi_BHBH(N, para)
    @unpack α,β,D,C,a,b,K,τ = para
    g = calc_g(para)
    #equilibrium point
    return 1-(1/(1+α+β*N))-(D*g)/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*g*N)
end

function findequil_BHBH(para)
    if BHBH_existence_check(para) > 0
        return 0.0 #no equilibrium point exists
    else
    return find_zero(N -> Nequi_BHBH(N, para), 0.5)
    end
end

function NdataBHBH_C(τrange,Cval, defaultpara)
    data = zeros(length(τrange))
    for i in eachindex(τrange)
        local_par = deepcopy(defaultpara)
        local_par.τ = τrange[i]
        local_par.C = Cval
        data[i] = findequil_BHBH(local_par)
    end
    return data
end

function NdataBHBH_β(τrange,βval, defaultpara)
    data = zeros(length(τrange))
    for i in eachindex(τrange)
        local_par = deepcopy(defaultpara)
        local_par.τ = τrange[i]
        local_par.β = βval
        data[i] = findequil_BHBH(local_par)
    end
    return data
end

let 
    τrange = 1:1:15
    Cval01data=NdataBHBH_C(τrange, 0.1, BevHoltPar(α=0.1,β=0.3,D=0.1,a=5.0,b=200.0,K=1.0))
    Cval05data=NdataBHBH_C(τrange, 0.5, BevHoltPar(α=0.1,β=0.3,D=0.1,a=5.0,b=200.0,K=1.0))
    Cval10data=NdataBHBH_C(τrange, 1.0, BevHoltPar(α=0.1,β=0.3,D=0.1,a=5.0,b=200.0,K=1.0))
    βval01data=NdataBHBH_β(τrange, 0.1, BevHoltPar(α=0.1,β=0.3,D=0.1,a=5.0,b=200.0,K=1.0))
    βval05data=NdataBHBH_β(τrange, 0.5, BevHoltPar(α=0.1,β=0.3,D=0.1,a=5.0,b=200.0,K=1.0))
    βval10data=NdataBHBH_β(τrange, 1.0, BevHoltPar(α=0.1,β=0.3,D=0.1,a=5.0,b=200.0,K=1.0))
    p1=scatter(τrange .+ randn(length(τrange)) .* 0.05, Cval01data, color="#FDE725FF",markersize=6, label="C=0.1")
    plot!(collect(1:1:15), Cval01data, label="",linestyle=:dash, color="#FDE725FF")
    scatter!(τrange .+ randn(length(τrange)) .* 0.05, Cval05data, color="#238A8DFF", markersize=6,label="C=0.5")
    plot!(collect(1:1:15), Cval05data, label="",linestyle=:dash, color="#238A8DFF")
    scatter!(τrange .+ randn(length(τrange)) .* 0.05, Cval10data, color="#440154FF", markersize=6,label="C=1.0")
    plot!(collect(1:1:15), Cval10data, label="",linestyle=:dash, color="#440154FF")
    xlabel!("τ")
    ylabel!("N*(τ)")
    p2=scatter(τrange .+ randn(length(τrange)) .* 0.05, βval01data, color="#FDE725FF",markersize=6, label="β=0.1")
    plot!(collect(1:1:15), βval01data, label="",linestyle=:dash, color="#FDE725FF")
    scatter!(τrange .+ randn(length(τrange)) .* 0.05, βval05data, color="#238A8DFF", markersize=6,label="β=0.5")
    plot!(collect(1:1:15), βval05data, label="",linestyle=:dash, color="#238A8DFF")
    scatter!(τrange .+ randn(length(τrange)) .* 0.05, βval10data, color="#440154FF", markersize=6,label="β=1.0")
    plot!(collect(1:1:15), βval10data, label="",linestyle=:dash, color="#440154FF")
    xlabel!("τ")
    ylabel!("N*(τ)")
    plot(p1, p2, layout=(2,1), size=(600,600), legend=:topright)
    savefig(joinpath(abpath(), "figs/BevHoltBevHolt_tauequi.pdf"))
end



let 
    time = 500
    timeseries = model_recursion(10.0,time,BevHoltPar(τ=1,α=0.1,β=0.3,C=0.3,D=0.1,a=15.0,b=200.0,K=1.0), BevertonHolt_modelII)
    println(timeseries[end])
    println(findequil_BHBH(BevHoltPar(τ=1,α=0.1,β=0.3,C=0.3,D=0.1,a=15.0,b=200.0,K=1.0)))
    plot(0:1:time,timeseries)
    xlims!(0,50)
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

