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





#Cohort dependent survival of immature individuals (immature individuals exposed to density effects within cohort - but not density effect with mature individuals)



let 
    time = 500
    timeseries = model_recursion(10.0,time,BevHoltPar(τ=13,α=0.3,β=0.3,C=0.3,D=0.3,a=15.0,b=200.0,K=1.0), BevertonHolt_modelII)
    println(timeseries[end])
    # println(findequil_BHBH(BevHoltPar(τ=1,α=0.1,β=0.3,C=0.3,D=0.1,a=15.0,b=200.0,K=1.0)))
    plot(0:1:time,timeseries)
    # xlims!(0,50)
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

