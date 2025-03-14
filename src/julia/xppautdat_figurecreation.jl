include("packages.jl")
include("TradeOffs_CommonCode.jl")

"""
    cleanxppautdat(file_path::String) -> DataFrame

Reads and processes data from a file, returning a cleaned DataFrame.

# Arguments
- `file_path::String`: Path to the data file.

# Returns
- `DataFrame`: Cleaned DataFrame with columns `:a`, `:Lowp`, and `:PointTypeName`.
"""
function cleanxppautdat(file_path)
    datalm = readdlm(file_path)
    data = DataFrame(datalm, [:a, :Lowp, :Highp, :PointType1, :LineNum, :PointType2])
    @select!(data, :a, :Lowp, :PointType1, :PointType2)
    @transform!(data, :PointType1 = Int.(:PointType1), :PointType2 = Int.(:PointType2))
    @transform!(data, :PointType = string.(:PointType1) .* string.(:PointType2))
    @transform!(data, :PointTypeName = ifelse.(:PointType .== "25", "BP", ifelse.(:PointType .== "23", "HP", "Other")))
    @select!(data, :a, :Lowp, :PointTypeName)
    data = sort(data, :Lowp)
    return data
end

let 
    alow=alowerconstraint(RickerPar(p=0.3,τ=2.0, α=0.1, β=0.3, b=200, K=1.0))
    arange=alow:0.1:15.0
    fillbottom = zeros(length(arange))
    data=cleanxppautdat("src/xppaut/RickerConstanttau2.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh=[phigherconstraint(RickerPar(a=aval,τ=2.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    plot(bp_data.a, bp_data.Lowp, label="Transcritical",lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker",lw=2)
    plot!(arange, fillbottom, fillrange = phigh, fillalpha = 0.2, c = 1, label = "Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(9.8, 14.5)
    ylims!(0.0, 1.0)
    title!("τ = 2.0")
end

let 
    alow=alowerconstraint(RickerPar(p=0.3,τ=3.0, α=0.1, β=0.3, b=200, K=1.0))
    arange=alow:0.1:15.0
    fillbottom = zeros(length(arange))
    data=cleanxppautdat("src/xppaut/RickerConstanttau3woPD.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh=[phigherconstraint(RickerPar(a=aval,τ=3.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    pddata=perioddoublecurve(arange, 3)
    # return pddata
    plot(bp_data.a, bp_data.Lowp, label="Transcritical",lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker",lw=2)
    plot!(pddata[:,1], pddata[:,2], label="Period doubling",lw=2)
    plot!(arange, fillbottom, fillrange = phigh, fillalpha = 0.2, c = 1, label = "Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(0.0, 15.0)
    ylims!(0.0, 1.0)
    title!("τ = 3.0")
end

let 
    alow=alowerconstraint(RickerPar(p=0.3,τ=4.0, α=0.1, β=0.3, b=200, K=1.0))
    arange=alow:0.1:15.0
    fillbottom = zeros(length(arange))
    data=cleanxppautdat("src/xppaut/RickerConstanttau4.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh=[phigherconstraint(RickerPar(a=aval,τ=4.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    plot(bp_data.a, bp_data.Lowp, label="Transcritical",lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker",lw=2)
    plot!(arange, fillbottom, fillrange = phigh, fillalpha = 0.2, c = 1, label = "Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(0.0, 15.0)
    ylims!(0.0, 1.0)
    title!("τ = 4.0")
end

let 
    alow=alowerconstraint(RickerPar(p=0.3,τ=5.0, α=0.1, β=0.3, b=200, K=1.0))
    arange=alow:0.1:22.0
    fillbottom = zeros(length(arange))
    data=cleanxppautdat("src/xppaut/RickerConstanttau5.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh=[phigherconstraint(RickerPar(a=aval,τ=5.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    plot(bp_data.a, bp_data.Lowp, label="Transcritical",lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker",lw=2)
    plot!(arange, fillbottom, fillrange = phigh, fillalpha = 0.2, c = 1, label = "Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(0.0, 22.0)
    ylims!(0.0, 1.0)
    title!("τ = 5.0")
end

let 
    datatau2 = cleanxppautdat("src/xppaut/RickerConstanttau2.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Lowp)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Lowp)
    datatau3 = cleanxppautdat("src/xppaut/RickerConstanttau3woPD.dat")    
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Lowp)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Lowp)
    datatau4 = cleanxppautdat("src/xppaut/RickerConstanttau4.dat")    
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Lowp)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Lowp)
    datatau5 = cleanxppautdat("src/xppaut/RickerConstanttau5.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Lowp)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Lowp)
    plot(bp_datatau2.a, bp_datatau2.Lowp, lw=2, linestyle=:solid, color=:blue,label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Lowp, lw=2, linestyle=:dash, color=:blue,label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Lowp, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Lowp, lw=2, linestyle=:dash, color=:purple,label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Lowp, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Lowp, lw=2, linestyle=:dash, color=:red,label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Lowp, lw=2, linestyle=:solid,color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Lowp, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    
    ylims!(0.0, 1.0)
    xlabel!("a")
    ylabel!("p")
end


#Tau create and then kill oscillations
let
    par=RickerPar(a=13.0,p=0.57205,τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50;upperval=5))
    p2=scatter(RCorbitdata[1], log.(RCorbitdata[2]),color=:black, label="")
    xlabel!("τ")
    ylabel!("log(N)")
    xlims!(0.0,5.0)
    # savefig(joinpath(abpath(), "figs/tauorbitdiagram_RickerRicker.pdf"))
end


#Tau allows equilibrium to exist and then kills stable point
datatau5 = cleanxppautdat("src/xppaut/RickerConstanttau5.dat")    
@subset(datatau5, :PointTypeName .== "BP")
@subset(datatau5, :a .<10.1)

let
    par=RickerPar(a=10.0138,p=0.464,τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50;upperval=5))
    p2=scatter(RCorbitdata[1], RCorbitdata[2],color=:black, label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0,5.0)
    # savefig(joinpath(abpath(), "figs/tauorbitdiagram_RickerRicker.pdf"))
end


