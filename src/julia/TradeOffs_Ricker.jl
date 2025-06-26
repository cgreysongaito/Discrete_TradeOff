include("packages.jl")
include("TradeOffs_CommonCode.jl")

#TODO LIST
#make function that calculates lower and upperbounds of a for all tau
#make function that calculates lower and upper bounds of p for all tau

##Ricker Model

#Density independent survival of immature individuals
#Check that beta can remove the bell shape in Ricker constant
let #figure of RickerConstant tau=3 showing two period doubling to something?   par=RickerPar(p=0.4, a=0.5,b=200.0,α=1.5,β=0.6)
    println(alowerconstraint(par))
    println(ahigherconstraint(par))
    data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; upperval=8))
    scatter(data[1], data[2],color=:black)
    xlabel!("τ")
    ylabel!("N")
end

let 
    tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=3,p=0.7, α=0.1, β=0.3, b=200, K=1.0), 500))
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("a")
    ylabel!("N")
    xlims!(7.3,7.8)
    plot(p3, layout=(1,1), size = (500,400), legend=false, guidefontsize=10, ms=2)
    savefig(joinpath(abpath(), "figs/aorbitdiagram_Rickerconstanttau3.pdf"))
end


let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=2), 50))
    # tau3data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=3), 50))
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
    savefig(joinpath(abpath(), "figs/aorbitdiagram_Rickerconstant.pdf"))
end

let 
    time = 10000
    finalts=9900
    timeseries = model_recursion(0.1,time,RickerPar(τ=1,p=0.6, a=29.839,α=0.1,β=0.3,K=1.0,b=200), RickerConstant_model)
    # return timeseries[finalts:end]
    p1=scatter(timeseries[finalts-1:end-1], timeseries[finalts:end])
    xlabel!("N(t-3)")
    ylabel!("N(t)")
    p2=scatter(0.0:1.0:101, timeseries[finalts:end])
    xlabel!("t")
    ylabel!("N(t)")
    plot(p1,p2, layout=(2,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
    # savefig(joinpath(abpath(), "figs/timeemedding_a108_Rickerconstant.pdf"))
end

let 
    time = 10000
    finalts=9900
    timeseries = model_recursion(0.1,time,RickerPar(τ=1,p=0.57205, a=13.0,α=0.1,β=0.3,K=1.0,b=200), RickerConstant_wofec_model, optτ=2)
    return timeseries
    # p1=scatter(timeseries[finalts-1:end-1], timeseries[finalts:end])
    # xlabel!("N(t-3)")
    # ylabel!("N(t)")
    # p2=scatter(0.0:1.0:101, timeseries[finalts:end])
    # xlabel!("t")
    # ylabel!("N(t)")
    # plot(p1,p2, layout=(2,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
    # # savefig(joinpath(abpath(), "figs/timeemedding_a108_Rickerconstant.pdf"))
end

let 
    time = 100000
    finalts=99000
    timeseries = model_recursion(0.1,time,RickerPar(τ=1,p=0.6, a=29.844,α=0.1,β=0.3,K=1.0,b=200), RickerConstant_model)
    # return timeseries[finalts:end]
    p1=scatter(timeseries[finalts-1:end-1], timeseries[finalts:end])
    xlabel!("N(t-3)")
    ylabel!("N(t)")
    p2=scatter(0.0:1.0:101, timeseries[finalts:end])
    xlabel!("t")
    ylabel!("N(t)")
    plot(p1,p2, layout=(2,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
    # savefig(joinpath(abpath(), "figs/timeemedding_a108_Rickerconstant.pdf"))
end

ahigherconstraint(RickerPar(τ=1,p=0.6, a=7.611502,α=0.1, β=0.3, b=200, K=1.0))

let 
    time = 1000000
    finalts=999500
    timeseries = model_recursion(19.51518,time,RickerPar(τ=3,p=0.7, a=7.611502,α=0.1, β=0.3, b=200, K=1.0), RickerConstant_model)
    timeseries2 = model_recursion(timeseries[end],time,RickerPar(τ=3,p=0.7, a=7.611502,α=0.1, β=0.3, b=200, K=1.0), RickerConstant_model)
    timeseries3 = model_recursion(timeseries2[end],time,RickerPar(τ=3,p=0.7, a=7.611502,α=0.1, β=0.3, b=200, K=1.0), RickerConstant_model)

    # p1=scatter(timeseries[finalts-3:end-3], timeseries[finalts:end])
    # xlabel!("N(t-3)")
    # ylabel!("N(t)")
    p2=scatter(0.0:1.0:501, timeseries3[finalts:end])
    xlabel!("t")
    ylabel!("N(t)")
    plot(p2, layout=(2,1), size = (500,400), legend=false, guidefontsize=10, ms=2)
    # savefig(joinpath(abpath(), "figs/timeemedding_a108_Rickerconstant.pdf"))
end


let 
    time = 5000
    finalts=4500
    timeseries = model_recursion(0.1,time,RickerPar(τ=3,p=0.6, a=11.0), RickerConstant_model)
    p1=scatter(timeseries[finalts-3:end-3], timeseries[finalts:end])
    xlabel!("N(t-3)")
    ylabel!("N(t)")
    p2=scatter(0.0:1.0:501, timeseries[finalts:end])
    xlabel!("t")
    ylabel!("N(t)")
    p3=scatter(0.0:1.0:51, timeseries[4950:end])
    xlabel!("t")
    ylabel!("N(t)")
    plot(p1,p2,p3, layout=(3,1), size = (600,900), legend=false, guidefontsize=10, ms=2)
    # savefig(joinpath(abpath(), "figs/timeemedding_a11_Rickerconstant.pdf"))
end

let 
    time = 5000
    finalts=4800
    timeseries = model_recursion(0.1,time,RickerPar(τ=3,p=0.6, a=11.0), RickerConstant_model)
    plot(timeseries[finalts-3:end-3], timeseries[finalts:end])
    xlabel!("N(t-3)")
    ylabel!("N(t)")
    # return timeseries[4502:5001]
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
    # savefig(joinpath(abpath(), "figs/aorbitdiagram_Rickerconstant_higherp.pdf"))
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
    savefig(joinpath(abpath(), "figs/porbitdiagram_Rickerconstant.pdf"))
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
    savefig(joinpath(abpath(), "figs/alphaorbitdiagram_Rickerconstant.pdf"))
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
    savefig(joinpath(abpath(), "figs/betaaorbitdiagram_Rickerconstant.pdf"))
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

#Using stability of first iterate to find period doubling (first Attempt) - 2nd attempt is to use the 2nd iterate


function firstiteratejacobian(τval, para)
    local_par = deepcopy(para)
    local_par.τ = τval
    @unpack a, b, K, p, τ = local_par
    matrix = zeros(Float64,τval+1,τval+1)
    matrix[1,1] = (1-calc_m(local_par))*(1+log(1-calc_m(local_par))+local_par.α)
    matrix[1,τval+1] = calc_m(local_par)
    for i in 2:τval+1
        matrix[i,i-1] = 1
    end
    return matrix
end #TODO update this to 

function firstiteratejacobian(τval, para)
    local_par = deepcopy(para)
    local_par.τ = τval
    @unpack a, b, K, p, τ, α,β = local_par
    matrix = zeros(Float64,τval+1,τval+1)
    matrix[1,1] = exp(-α-β*x)*(1-β*x)
    matrix[1,τval+1] = calc_m(local_par)
    for i in 2:τval+1
        matrix[i,i-1] = 1
    end
    return matrix
end


function eigenvals(pval, τval)
    arange=3.0:0.0001:9.0
    eig1=[]
    eig2=[]
    eig3=[]
    eig4=[]
    maxeig=[]
    endai=0
    for ai in eachindex(arange)
        local_par=RickerPar(τ=τval, a=arange[ai],α=0.1,β=0.3,b=200,K=1.0,p=pval)
        if calc_m(local_par)>=1.0
            endai=ai-1
            break
        else
            matrix = firstiteratejacobian(τval, local_par)
            eigs = real(eigen(matrix).values)
            maxeigval=maximum(abs.(eigs))
            push!(eig1, eigs[1])
            push!(eig2, eigs[2])
            push!(eig3, eigs[3])
            push!(eig4, eigs[4])
            push!(maxeig, maxeigval)
        end
    end
    return [arange[1:endai],eig1,eig2,eig3,eig4, maxeig]
end

function findflip(eigdata)
    data=[]
    for i in 2:length(eigdata)-1
        for j in 1:length(eigdata[i])
            if isapprox(eigdata[i][j], 1.0,atol=1e-5) || isapprox(eigdata[i][j], -1.0,atol=1e-5)
                push!(data, eigdata[1][j])
            end
        end
    end
    return data
end
findflip(eigenvals(0.7, 3))

test =eigenvals(0.7, 3)
isapprox(test[2][3],-0.25,atol=1e-3)
test[2][end]
test[1][end]
let 
    eigs = eigenvals(0.7, 3)
    plot(eigs[1],eigs[2],label="eig1")
    plot!(eigs[1],eigs[3],label="eig2")
    plot!(eigs[1],eigs[4],label="eig3")
    plot!(eigs[1],eigs[5],label="eig4")
    hline!([1.0], linestyle=:dash, color=:black, linewidth=1.5, label="")
    hline!([-1.0], linestyle=:dash, color=:black, linewidth=1.5, label="")
    xlabel!("a")
    ylabel!("λ")
    title!("p=0.7, τ=3")
end

function seconditeratejacobian(τval, para)
    local_par = deepcopy(para)
    local_par.τ = τval
    m=calc_m(local_par)
    @unpack α = local_par
    matrix = zeros(Float64,τval+1,τval+1)
    matrix[1,1] = ((1-m)^2)*(1+log(1-m)+α)^2
    matrix[1,3] = m
    matrix[1,τval+1] = m*(1-m)*(1+log(1-m)+α)
    matrix[2,1] = (1-m)*(1+log(1-m)+α)
    matrix[2,4] = m
    for i in 3:τval+1
        matrix[i,i-1] = 1
    end
    return matrix
end

test=seconditeratejacobian(3, RickerPar(a=10.0,α=0.1,β=0.3,b=200,K=1.0,p=0.6))
eigen(test).values

function eigenvals2nditerate(aval, τval)
    prange=0.1:0.00001:1.0
    eig1=[]
    eig2=[]
    eig3=[]
    eig4=[]
    maxeig=[]
    endpi=0
    for pi in eachindex(prange)
        local_par=RickerPar(τ=τval, a=aval,α=0.1,β=0.3,b=200,K=1.0,p=prange[pi])
        if calc_m(local_par)>=1.0
            endpi=pi-1
            break
        else
            matrix = seconditeratejacobian(τval, local_par)
            eigs = real(eigen(matrix).values)
            maxeigval=maximum(abs.(eigs))
            push!(eig1, eigs[1])
            push!(eig2, eigs[2])
            push!(eig3, eigs[3])
            push!(eig4, eigs[4])
            push!(maxeig, maxeigval)
        end
    end
    return [prange[1:endpi],eig1,eig2,eig3,eig4, maxeig]
end

let 
    eigs = eigenvals2nditerate(10.0, 3)
    plot(eigs[1],eigs[2],label="eig1")
    plot!(eigs[1],eigs[3],label="eig2")
    plot!(eigs[1],eigs[4],label="eig3")
    plot!(eigs[1],eigs[5],label="eig4")
    hline!([1.0], linestyle=:dash, color=:black, linewidth=1.5, label="")
    hline!([-1.0], linestyle=:dash, color=:black, linewidth=1.5, label="")
    xlabel!("p")
    ylabel!("λ")
    title!("a=10.0, τ=3")
end

function findperioddouble(aval, τval)
    if aval<5.0
        prange=0.7:0.00001:1.0
    else
        prange=0.52:0.00001:1.0
    end
    for pi in eachindex(prange)
        matrix = firstiteratejacobian(τval, RickerPar(a=aval,α=0.1,β=0.3,b=200,K=1.0,p=prange[pi]))
        domeig=maximum(abs.(real(eigen(matrix).values)))
        if domeig>1.0 && !isapprox(domeig, 1.0)
            return prange[pi]
        end
    end
end

function perioddoublecurve(arange, τval)
    pdata = Vector{Union{Float64, Nothing}}(undef, length(arange))
    adata = Vector{Union{Float64, Nothing}}(undef, length(arange))
    @threads for ai in eachindex(arange)
        pdata[ai] = findperioddouble(arange[ai], τval)
        adata[ai] = arange[ai]
    end
    combined_data = hcat(adata, pdata)
    filtered_data = combined_data[.!isnothing.(pdata), :]
    return filtered_data
end

perioddoublecurve(4.8:0.1:15.0, 3)

#symbolic differentiation
using SymPy
@syms x y z
@syms α β m
f(x,y,z)= exp(-α-β*(exp(-α-β*x)*x+m*z))*(exp(-α-β*x)*x+m*z)+m*y
diff(f(x,y,z),x)
diff(f(x,y,z),y)
diff(f(x,y,z),z)

g(x)=exp(-α-β*(exp(-α-β*x)*x+m*((exp(-α-β*x)*x)/(1-m))))*(exp(-α-β*x)*x+m*((exp(-α-β*x)*x)/(1-m)))+m*x-x

solve(exp(-α-β*(exp(-α-β*x)*x+m*((exp(-α-β*x)*x)/(1-m))))*(exp(-α-β*x)*x+m*((exp(-α-β*x)*x)/(1-m)))+m*x-x)

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
    savefig(joinpath(abpath(), "figs/aorbitdiagram_RickerBevertonDefault.pdf"))
end

let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=2,α=0.0005, D=0.0005), 50; upperval=60))
    tau3data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=3,α=0.0005, D=0.0005), 50; upperval=60))
    tau4data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=4,α=0.0005, D=0.0005), 50; upperval=60))
    tau5data=flattenorbitdata(orbitdiagrams(RickerBeverton_model, "a", RickerPar(τ=5,α=0.0005, D=0.0005), 50; upperval=60))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,60.0)
    p3=scatter(tau3data[1], tau3data[2],color=:black)
    title!("τ=3")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,60.0)
    p4=scatter(tau4data[1], tau4data[2],color=:black)
    title!("τ=4")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,60.0)
    p5=scatter(tau5data[1], tau5data[2],color=:black)
    title!("τ=5")
    xlabel!("a")
    ylabel!("N")
    xlims!(0.0,60.0)
    plot(p2,p3,p4,p5, layout=(4,1), size = (500,700), legend=false, guidefontsize=10, ms=2)
    savefig(joinpath(abpath(), "figs/aorbitdiagram_RickerBevertonlowalphalowD.pdf"))
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
    time = 5000
    τval=15
    timeseries = model_Leslierecursion(τval, time, RickerPar(τ=τval, a=100.0,α=1.5,β=0.1,D=0.9,C=0.15), 0.1, LeslieMatrix)
    # plot(0:1:time,first_elements(timeseries))
    return first_elements(timeseries)[1:1200]
end


calc_g(RickerPar(τ=1, a=25.0,α=1.5,β=2.0,D=0.5,C=0.15))

function LeslieMatrixOrbitDiagram(τrange, time, finalts, para, init, lesliematrix)
    data = Vector{Vector{Float64}}(undef, length(τrange))
    @threads for τi in eachindex(τrange)
        timeseries = model_Leslierecursion(τrange[τi], time, para, init, lesliematrix)
        data[τi] = first_elements(timeseries)[end-finalts:end]
    end
    return [τrange,data]
end

let 
    timeseries = model_recursion(0.1,5000,RickerPar(τ=0, a=100.0,α=1.5,β=2.0,D=0.5,C=0.15), RickerLeslie_τ0_model)
    # return timeseries[end-100:end]
    scatter(0.0:1.0:100.0, timeseries[end-100:end])
end



let
    τrange = 1:1:15
    # RIorbitdata = model_τorbit(τrange, RickerConstant_model, RickerPar(a=100, p=0.6), 50)
    # RIIorbitdata = model_τorbit(τrange, RickerBeverton_model, RickerPar(a=100,D=0.5,C=0.15), 50)
    par=RickerPar(τ=0, a=100.0,α=1.5,β=0.1,D=0.9,C=0.15)
    RLorbitdata = flattenorbitdata(LeslieMatrixOrbitDiagram(τrange, 50000, 100, par, 0.1, LeslieMatrix))
    
    return RL
    # RLτ0data = model_recursion(0.1,5000,par, RickerLeslie_τ0_model)
    # RLτ0data_trans=RLτ0data[end-100:end]
    # p2=scatter(RLorbitdata[1], RLorbitdata[2],color=:black, label="")
    # scatter!(zeros(length(RLτ0data_trans)), RLτ0data_trans, color=:black, label="")
    # xlabel!("τ")
    # ylabel!("N")
    # xlims!(0.0,15.0)
    # savefig(joinpath(abpath(), "figs/tauorbitdiagram_RickerRicker.pdf"))
end

calc_g(RickerPar(τ=0, a=74,α=1.5,β=2.0,D=0.5,C=0.15))
alowerconstraint(RickerPar(τ=0, a=10.0,α=1.5,β=2.0,D=0.5,C=0.15))
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


#Xppaut helper
tau2par=RickerPar(τ=2.0, a=5.0,α=0.1,β=0.3,b=200,K=1.0,p=0.6)
alowerconstraint(tau2par)
ahigherconstraint(tau2par)
calc_m(tau2par)

tau3par=RickerPar(τ=3.0, a=5.0,α=0.1,β=0.3,b=200,K=1.0)
alowerconstraint(tau3par)

tau4par=RickerPar(τ=4.0, a=5.0,α=0.1,β=0.3,b=200,K=1.0,p=0.6)
alowerconstraint(tau4par)
ahigherconstraint(tau4par)

tau5par=RickerPar(τ=5.0, a=5.0,α=0.1,β=0.3,b=200,K=1.0,p=0.6)
alowerconstraint(tau5par)
ahigherconstraint(tau5par)

let 
    tau2data=flattenorbitdata(orbitdiagrams(RickerConstant_model, "a", RickerPar(τ=2,p=0.6,α=0.1,β=0.3,b=200,K=1.0), 50))
    p2=scatter(tau2data[1], tau2data[2],color=:black)
    title!("τ=2")
    xlabel!("a")
    ylabel!("N")
    xlims!(10.0,14.5)
    ylims!(0.0,15.0)
end