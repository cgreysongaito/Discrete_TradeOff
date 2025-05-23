include("packages.jl")
default(grid=false, linewidth=3, tickfontsize=12, legendfontsize=10, guidefontsize=15)
include("TradeOffs_CommonCode.jl")

"""
    cleanxppautdat(file_path::String) -> DataFrame

Reads and processes data from a file, returning a cleaned DataFrame.

# Arguments
- `file_path::String`: Path to the data file.

# Returns
- `DataFrame`: Cleaned DataFrame with columns `:a`, `:Lowp`, and `:PointTypeName`.
"""
function cleanxppautdat_onepar(file_path)
    datalm = readdlm(file_path)
    data = DataFrame(datalm, [:a, :N1, :N2, :PointType1, :LineNum, :PointType2])
    @select!(data, :a, :N1, :PointType1, :PointType2, :LineNum)
    @transform!(data, :PointType1 = Int.(:PointType1), :PointType2 = Int.(:PointType2))
    @transform!(data, :PointType = string.(:PointType1) .* string.(:PointType2))
    @transform!(data, :PointTypeName = ifelse.(:PointType .== "10", "Stable", ifelse.(:PointType .== "20", "Unstable", "Other")))
    @select!(data, :a, :N1, :PointTypeName, :LineNum)
    # data = sort(data, :N1)
    return data
end

function cleanxppautdat_twopar(file_path)
    datalm = readdlm(file_path)
    data = DataFrame(datalm, [:a, :Low2ndpar, :High2ndpar, :PointType1, :LineNum, :PointType2])
    @select!(data, :a, :Low2ndpar, :PointType1, :PointType2, :LineNum)
    @transform!(data, :PointType1 = Int.(:PointType1), :PointType2 = Int.(:PointType2))
    @transform!(data, :PointType = string.(:PointType1) .* string.(:PointType2))
    @transform!(data, :PointTypeName = ifelse.(:PointType .== "25", "BP", ifelse.(:PointType .== "23", "HP", "Other")))
    @select!(data, :a, :Low2ndpar, :PointTypeName, :LineNum)
    data = sort(data, :Low2ndpar)
    return data
end

#RickerConstant
let #Even tau
    datatau2 = cleanxppautdat_onepar("src/xppaut/RickerConstanttau2_a.dat")
    tau2lowerbound = alowerconstraint(RickerPar(p=0.6, τ=2.0, α=0.1, β=0.3, b=200, K=1.0))
    datatau2_filtered = sort(@subset(datatau2, :a .> tau2lowerbound), :a) #.|| :N1 .>= 0
    datatau2s = @subset(datatau2_filtered, :PointTypeName .== "Stable" .&& :N1 .> 0.00)
    datatau2ua = @subset(datatau2_filtered, :PointTypeName .== "Unstable" .&& :N1 .> 6.00)
    datatau2ub = @subset(datatau2_filtered, :PointTypeName .== "Unstable" .&& :N1 .< 1.00)
    datatau4 = cleanxppautdat_onepar("src/xppaut/RickerConstanttau4_a.dat")
    tau4lowerbound = alowerconstraint(RickerPar(p=0.6, τ=4.0, α=0.1, β=0.3, b=200, K=1.0))
    datatau4_filtered = sort(@subset(datatau4, :a .> tau4lowerbound), :a) #.|| :N1 .>= 0
    datatau4s = @subset(datatau4_filtered, :PointTypeName .== "Stable" .&& :N1 .> 0.00)
    datatau4ua = @subset(datatau4_filtered, :PointTypeName .== "Unstable" .&& :N1 .> 6.00)
    datatau4ub = @subset(datatau4_filtered, :PointTypeName .== "Unstable" .&& :N1 .< 1.00)
    datatau6 = cleanxppautdat_onepar("src/xppaut/RickerConstanttau6_a.dat")
    tau6lowerbound = alowerconstraint(RickerPar(p=0.6, τ=6.0, α=0.1, β=0.3, b=200, K=1.0))
    datatau6_filtered = sort(@subset(datatau6, :a .> tau4lowerbound), :a) #.|| :N1 .>= 0
    datatau6s = unique(@subset(datatau6_filtered, :PointTypeName .== "Stable" .&& :N1 .> 0.00), :a)
    datatau6ua = @subset(datatau6_filtered, :PointTypeName .== "Unstable" .&& :N1 .> 6.00)
    datatau6ub = @subset(datatau6_filtered, :PointTypeName .== "Unstable" .&& :N1 .< 1.00)
    plot(datatau2s.a, datatau2s.N1, color=:black, label="τ=2.0")
    plot!(datatau2ua.a, datatau2ua.N1, color=:black, linestyle=:dash, label="")
    plot!(datatau4s.a, datatau4s.N1, color=:blue, label="τ=4.0")
    plot!(datatau4ua.a, datatau4ua.N1, color=:blue, linestyle=:dash, label="")
    plot!(datatau6s.a, datatau6s.N1, color=:red, label="τ=6.0")
    plot!(datatau6ua.a, datatau6ua.N1, color=:red, linestyle=:dash, label="")
    plot!([-1], [0], linestyle=:solid, color=:black, label="Stable")
    plot!([-1], [0], linestyle=:dash, color=:black, label="Unstable")
    xlabel!("a")
    ylabel!("N*")
    ylims!(-0.5, 15.0)
    xlims!(0.0, 38.0)
    savefig(joinpath(abpath(), "figs/RickerConstanttaueven_a.pdf"))
end

function branchsplitter2(subsetteddata)
    series = []
    current_series = DataFrame()
    n = nrow(subsetteddata)
    i = 1
    while i <= n
        push!(current_series, subsetteddata[i, :])
        if i < n && abs(subsetteddata.N1[i+1] - subsetteddata.N1[i]) > 0.5
            push!(series, current_series)
            current_series = DataFrame()
        end
        i += 1
    end
    if nrow(current_series) > 0
        push!(series, current_series)
    end
    return series
end

let #tau=3 (odd)
    datatau3 = cleanxppautdat_onepar("src/xppaut/RickerConstanttau3_a.dat")
    tau3lowerbound = alowerconstraint(RickerPar(p=0.6, τ=3.0, α=0.1, β=0.3, b=200, K=1.0))
    tau3upperbound = ahigherconstraint(RickerPar(p=0.6, τ=3.0, α=0.1, β=0.3, b=200, K=1.0))
    datatau3_filtered = @subset(datatau3, :a .> tau3lowerbound)
    datatau3s1 = @subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 1.0)
    datatau3s2p2upper = unique(@subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 2.0 .&& :N1 .> maximum(datatau3s1.N1)), :a)
    datatau3s2p2lower = unique(@subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 2.0 .&& :N1 .< maximum(datatau3s1.N1) .&& :N1 .> 0.1), :a)
    datatau3s2_0 = @subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 2.0 .&& :N1 .< 1.00)
    datatau3s3p4upperupper = @subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 3.0 .&& :N1 .> maximum(datatau3s2p2upper.N1))
    datatau3s3p4upperlower = @subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 3.0 .&& :N1 .< maximum(datatau3s2p2upper.N1) .&& :N1 .> 0.1)
    datatau3s3p4lowerupper = unique(@subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 4.0 .&& :N1 .< maximum(datatau3s2p2lower.N1)), :a)
    datatau3s3p4lowerlower = @subset(datatau3_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 4.0 .&& :N1 .< maximum(datatau3s2p2lower.N1) .&& :N1 .> 0.1)

    datatau3u1 = unique(@subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7 .&& :LineNum .== 1.0 .&& :N1 .> 0.00), :a)
    datatau3u2upper = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .> maximum(datatau3s1.N1) .&& :a .> maximum(datatau3s1.a) + 0.1)
    datatau3u2lower = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .< maximum(datatau3s1.N1) .&& :N1 .> 0 .&& :a .> maximum(datatau3s1.a) + 0.1)
    datatau3u3 = unique(@subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 3.0 .&& :N1 .> 0.00), :a)
    datatau3u4upper = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 4.0 .&& :N1 .> maximum(datatau3s3p4lowerlower.N1))
    datatau3u4lower = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 4.0 .&& :N1 .< maximum(datatau3s3p4lowerlower.N1) .&& :N1 .> 0.001)
    datatau3u0 = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .< 0.001)
    datatau3ul0 = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 1.0 .&& :N1 .< 0.001)

    plot(datatau3s1.a, datatau3s1.N1, color=:black, label="Stable")
    plot!(datatau3s2p2upper.a, datatau3s2p2upper.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3s2p2lower.a, datatau3s2p2lower.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3s2_0.a, datatau3s2_0.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3s3p4upperupper.a, datatau3s3p4upperupper.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3s3p4upperlower.a, datatau3s3p4upperlower.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3s3p4lowerlower.a, datatau3s3p4lowerlower.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3s3p4lowerupper.a, datatau3s3p4lowerupper.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau3u1.a, datatau3u1.N1, color=:black, linestyle=:dash, lw=1.5, label="Unstable")
    plot!(datatau3u2upper.a, datatau3u2upper.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u2lower.a, datatau3u2lower.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u3.a, datatau3u3.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u4upper.a, datatau3u4upper.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u4lower.a, datatau3u4lower.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u0.a, datatau3u0.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3ul0.a, datatau3ul0.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    ylabel!("N*")
    xlabel!("a")
    xlims!(0.0, tau3upperbound + 0.1)
    ylims!(-0.5, 30.0)
    println(maximum(datatau3s3p4lowerupper.a))
end

let #Time embedding for tau=3 when period 4
    time = 1000000
    finalts = 999000
    timeseries = model_recursion(0.1, time, RickerPar(τ=3, p=0.6, a=10.99, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    df = DataFrame(time = 0:time, N = timeseries)
    df.periodpoint = mod.(df.time, 4) .+ 1
    # timeseries1 = @subset(df, :periodpoint .== 1)
    # timeseries2 = @subset(df, :periodpoint .== 2)
    # timeseries3 = @subset(df, :periodpoint .== 3)
    # timeseries4 = @subset(df, :periodpoint .== 4)
    dfxaxis = @subset(df, :time .>= finalts-3)
    dfyaxis = @subset(df, :time .>= finalts)
    return
    scatter(dfxaxis.Nlt s-1:end-1], dfwotrans, color=dfwotrans.periodpoint)
    xlabel!("N(t-3)")
    ylabel!("N(t)")
    # savefig(joinpath(abpath(), "figs/timeemedding_a108_Rickerconstant.pdf"))
end

let #Time embedding for tau=3 when N-S?

end

let
    alow = alowerconstraint(RickerPar(p=0.3, τ=2.0, α=0.1, β=0.3, b=200, K=1.0))
    arange = alow:0.1:15.0
    fillbottom = zeros(length(arange))
    data = cleanxppautdat("src/xppaut/RickerConstanttau2.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh = [phigherconstraint(RickerPar(a=aval, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
    plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(9.8, 14.5)
    ylims!(0.0, 1.0)
    title!("τ = 2.0")
end

let
    alow = alowerconstraint(RickerPar(p=0.3, τ=3.0, α=0.1, β=0.3, b=200, K=1.0))
    arange = alow:0.1:15.0
    fillbottom = zeros(length(arange))
    data = cleanxppautdat("src/xppaut/RickerConstanttau3woPD.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh = [phigherconstraint(RickerPar(a=aval, τ=3.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    pddata = perioddoublecurve(4.3:0.1:15.0, 3)
    # return pddata
    plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
    plot!(pddata[:, 1], pddata[:, 2], label="Period doubling", lw=2)
    plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(0.0, 15.0)
    ylims!(0.0, 1.0)
    title!("τ = 3.0")
    savefig(joinpath(abpath(), "figs/RickerConstanttau3wPD.pdf"))
end

let
    alow = alowerconstraint(RickerPar(p=0.3, τ=4.0, α=0.1, β=0.3, b=200, K=1.0))
    arange = alow:0.1:15.0
    fillbottom = zeros(length(arange))
    data = cleanxppautdat("src/xppaut/RickerConstanttau4.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    phigh = [phigherconstraint(RickerPar(a=aval, τ=4.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
    plot!(pddata[:, 1], pddata[:, 2], label="Period doubling", lw=2)
    plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
    xlabel!("a")
    ylabel!("p")
    xlims!(0.0, 15.0)
    ylims!(0.0, 1.0)
    title!("τ = 4.0")
end

let
    alow = alowerconstraint(RickerPar(p=0.3, τ=5.0, α=0.1, β=0.3, b=200, K=1.0))
    arange = alow:0.1:22.0
    fillbottom = zeros(length(arange))
    data = cleanxppautdat("src/xppaut/RickerConstanttau5.dat")
    bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
    hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
    pddata = perioddoublecurve(5.6:0.1:22.0, 5)
    phigh = [phigherconstraint(RickerPar(a=aval, τ=5.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
    plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
    plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
    plot!(pddata[:, 1], pddata[:, 2], label="Period doubling", lw=2)
    plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
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
    plot(bp_datatau2.a, bp_datatau2.Lowp, lw=2, linestyle=:solid, color=:blue, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Lowp, lw=2, linestyle=:dash, color=:blue, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Lowp, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Lowp, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Lowp, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Lowp, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Lowp, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Lowp, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 22.0)
    ylims!(0.0, 1.0)
    xlabel!("a")
    ylabel!("p")
    savefig(joinpath(abpath(), "figs/apbifurcation_RickerConstant.pdf"))
end


#Tau create and then kill oscillations
let
    par = RickerPar(a=13.0, p=0.57205, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; upperval=5))
    p2 = scatter(RCorbitdata[1], log.(RCorbitdata[2]), color=:black, label="")
    xlabel!("τ")
    ylabel!("log(N)")
    xlims!(0.0, 5.0)
    # savefig(joinpath(abpath(), "figs/tauorbitdiagram_RickerRicker.pdf"))
end


#Tau allows equilibrium to exist and then kills stable point
datatau5 = cleanxppautdat("src/xppaut/RickerConstanttau5.dat")
@subset(datatau5, :PointTypeName .== "BP")
@subset(datatau5, :a .< 10.1)

let
    par = RickerPar(a=10.0138, p=0.464, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; upperval=5))
    p2 = scatter(RCorbitdata[1], RCorbitdata[2], color=:black, label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 5.0)
    # savefig(joinpath(abpath(), "figs/tauorbitdiagram_RickerRicker.pdf"))
end


#RickerRicker
#a and C two par bifurcation figure

let
    datatau1 = cleanxppautdat("src/xppaut/RickerRickertau1_C.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat("src/xppaut/RickerRickertau2_C.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/RickerRickertau3_C.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/RickerRickertau4_C.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/RickerRickertau5_C.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=1")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 50.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("C")
    savefig(joinpath(abpath(), "figs/aCbifurcation_RickerRicker.pdf"))
end

let
    datatau1 = cleanxppautdat("src/xppaut/RickerRickertau1_C.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat("src/xppaut/RickerRickertau2_C.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/RickerRickertau3_C.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/RickerRickertau4_C.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/RickerRickertau5_C.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=1")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 50.0)
    ylims!(0.0, 0.6)
    xlabel!("a")
    ylabel!("C")
    title!("Zoomed in")
    savefig(joinpath(abpath(), "figs/aCbifurcation_RickerRicker_weakcomp.pdf"))
end

let
    datatau1 = cleanxppautdat("src/xppaut/RickerRickertau1_beta.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat("src/xppaut/RickerRickertau2_beta.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/RickerRickertau3_beta.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/RickerRickertau4_beta.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/RickerRickertau5_beta.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=1")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 70.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("β")
    savefig(joinpath(abpath(), "figs/abetabifurcation_RickerRicker.pdf"))
end

let
    datatau1 = cleanxppautdat("src/xppaut/RickerRickertau1_beta.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat("src/xppaut/RickerRickertau2_beta.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/RickerRickertau3_beta.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/RickerRickertau4_beta.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/RickerRickertau5_beta.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=1")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 70.0)
    ylims!(0.0, 0.3)
    xlabel!("a")
    ylabel!("β")
    title!("Zoomed in")
    savefig(joinpath(abpath(), "figs/abetabifurcation_RickerRicker_weak.pdf"))
end

#BevertonHoltRicker
#a and C two par bifurcation figure

let
    datatau0 = cleanxppautdat("src/xppaut/BevertonHoltRickertau0_C.dat")
    bp_datatau0 = unique(@subset(datatau0, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau0 = sort(unique(@subset(datatau0, :PointTypeName .== "HP"), :Low2ndpar), :a)
    datatau1 = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_C.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = sort(unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar), :a)
    datatau2 = cleanxppautdat("src/xppaut/BevertonHoltRickertau2_C.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/BevertonHoltRickertau3_C.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/BevertonHoltRickertau4_C.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/BevertonHoltRickertau5_C.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau0.a, bp_datatau0.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=0")
    plot!(hp_datatau0.a, hp_datatau0.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=0")
    plot!(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:brown, label="BP τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="HP τ=1")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:orange, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:orange, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 50.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("C")
    savefig(joinpath(abpath(), "figs/aCbifurcation_BevertonHoltRicker.pdf"))
end

let
    datatau0 = cleanxppautdat("src/xppaut/BevertonHoltRickertau0_C.dat")
    bp_datatau0 = unique(@subset(datatau0, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau0 = sort(unique(@subset(datatau0, :PointTypeName .== "HP"), :Low2ndpar), :a)
    datatau1 = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_C.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = sort(unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar), :a)
    datatau2 = cleanxppautdat("src/xppaut/BevertonHoltRickertau2_C.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/BevertonHoltRickertau3_C.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/BevertonHoltRickertau4_C.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/BevertonHoltRickertau5_C.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau0.a, bp_datatau0.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=0")
    plot!(hp_datatau0.a, hp_datatau0.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=0")
    plot!(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:brown, label="BP τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="HP τ=1")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:orange, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:orange, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 200.0)
    ylims!(0.0, 0.15)
    xlabel!("a")
    ylabel!("C")
    title!("Zoomed in")
    savefig(joinpath(abpath(), "figs/aCbifurcation_BevertonHoltRicker_weakcomp.pdf"))
end

let
    datatau1 = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_beta.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = sort(unique(@subset(datatau1, :PointTypeName .== "HP" .&& :a .<= 50.0), :a), :a)
    CSV.write("hp_datatau1.csv", hp_datatau1)
    # return hp_datatau1
end


let
    datatau0 = cleanxppautdat("src/xppaut/BevertonHoltRickertau0_beta.dat")
    bp_datatau0 = unique(@subset(datatau0, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau0 = sort(unique(@subset(datatau0, :PointTypeName .== "HP"), :a), :a)
    datatau1_a = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_beta_adata.dat")
    datatau1_b = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_beta_bdata.dat")
    bp_datatau1 = unique(@subset(datatau1_a, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1a = unique(@subset(datatau1_a, :PointTypeName .== "HP"), :Low2ndpar)
    hp_datatau1b = sort(unique(@subset(datatau1_b, :PointTypeName .== "HP"), :a), :a)
    hp_datatau1 = vcat(hp_datatau1a, hp_datatau1b)
    datatau2 = cleanxppautdat("src/xppaut/BevertonHoltRickertau2_beta.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/BevertonHoltRickertau3_beta.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/BevertonHoltRickertau4_beta.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/BevertonHoltRickertau5_beta.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau0.a, bp_datatau0.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=0")
    plot!(hp_datatau0.a, hp_datatau0.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=0")
    plot!(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:brown, label="BP τ=1")
    plot!(hp_datatau1a.a, hp_datatau1a.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="HP τ=1")
    plot!(hp_datatau1b.a, hp_datatau1b.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:orange, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:orange, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0, 100.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("β")
    savefig(joinpath(abpath(), "figs/abetabifurcation_BevertonHoltRicker.pdf"))
end

let
    datatau0 = cleanxppautdat("src/xppaut/BevertonHoltRickertau0_beta.dat")
    bp_datatau0 = unique(@subset(datatau0, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau0 = sort(unique(@subset(datatau0, :PointTypeName .== "HP"), :a), :a)
    datatau1_a = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_beta_adata.dat")
    datatau1_b = cleanxppautdat("src/xppaut/BevertonHoltRickertau1_beta_bdata.dat")
    bp_datatau1 = unique(@subset(datatau1_a, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1a = unique(@subset(datatau1_a, :PointTypeName .== "HP"), :Low2ndpar)
    hp_datatau1b = sort(unique(@subset(datatau1_b, :PointTypeName .== "HP"), :a), :a)
    hp_datatau1 = vcat(hp_datatau1a, hp_datatau1b)
    datatau2 = cleanxppautdat("src/xppaut/BevertonHoltRickertau2_beta.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat("src/xppaut/BevertonHoltRickertau3_beta.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat("src/xppaut/BevertonHoltRickertau4_beta.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat("src/xppaut/BevertonHoltRickertau5_beta.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau0.a, bp_datatau0.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="BP τ=0")
    plot!(hp_datatau0.a, hp_datatau0.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="HP τ=0")
    plot!(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:brown, label="BP τ=1")
    plot!(hp_datatau1a.a, hp_datatau1a.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="HP τ=1")
    plot!(hp_datatau1b.a, hp_datatau1b.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:orange, label="BP τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:orange, label="HP τ=2")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="BP τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="HP τ=3")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="BP τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="HP τ=4")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="BP τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="HP τ=5")
    xlims!(0.0, 400.0)
    ylims!(0.0, 0.3)
    xlabel!("a")
    ylabel!("β")
    title!("Zoomed in")
    savefig(joinpath(abpath(), "figs/abetabifurcation_BevertonHoltRicker_weak.pdf"))
end




#BevertonHoltBevertonHolt model
let
    datatau0 = cleanxppautdat_onepar("src/xppaut/BevertonHoltBevertonHolttau0_a.dat")
    sta_data0 = @subset(datatau0, :PointTypeName .== "Stable")
    datatau1 = cleanxppautdat_onepar("src/xppaut/BevertonHoltBevertonHolttau1_a.dat")
    sta_data1 = @subset(datatau1, :PointTypeName .== "Stable")
    datatau2 = cleanxppautdat_onepar("src/xppaut/BevertonHoltBevertonHolttau2_a.dat")
    sta_data2 = @subset(datatau2, :PointTypeName .== "Stable")
    datatau3 = cleanxppautdat_onepar("src/xppaut/BevertonHoltBevertonHolttau3_a.dat")
    sta_data3 = @subset(datatau3, :PointTypeName .== "Stable")
    datatau4 = cleanxppautdat_onepar("src/xppaut/BevertonHoltBevertonHolttau4_a.dat")
    sta_data4 = @subset(datatau4, :PointTypeName .== "Stable")
    datatau5 = cleanxppautdat_onepar("src/xppaut/BevertonHoltBevertonHolttau5_a.dat")
    sta_data5 = @subset(datatau5, :PointTypeName .== "Stable")
    plot(sta_data0.a, sta_data0.N1, lw=2, linestyle=:solid, color=:blue, label="τ=0")
    plot!(sta_data1.a, sta_data1.N1, lw=2, linestyle=:solid, color=:brown, label="τ=1")
    plot!(sta_data2.a, sta_data2.N1, lw=2, linestyle=:solid, color=:orange, label="τ=2")
    plot!(sta_data3.a, sta_data3.N1, lw=2, linestyle=:solid, color=:purple, label="τ=3")
    plot!(sta_data4.a, sta_data4.N1, lw=2, linestyle=:solid, color=:red, label="τ=4")
    plot!(sta_data5.a, sta_data5.N1, lw=2, linestyle=:solid, color=:green, label="τ=5")
    xlims!(0.0, 500.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("N")
    title!("BevertonHolt Beverton Holt model")
    savefig(joinpath(abpath(), "figs/BevertonHoltBevertonHolt_a.pdf"))
end

#RickerBevertonHolt model
let
    datatau0 = cleanxppautdat_onepar("src/xppaut/RickerBevertonHolttau0_a.dat")
    sta_data0 = @subset(datatau0, :PointTypeName .== "Stable")
    datatau1 = cleanxppautdat_onepar("src/xppaut/RickerBevertonHolttau1_a.dat")
    sta_data1 = @subset(datatau1, :PointTypeName .== "Stable")
    datatau2 = cleanxppautdat_onepar("src/xppaut/RickerBevertonHolttau2_a.dat")
    sta_data2 = @subset(datatau2, :PointTypeName .== "Stable")
    datatau3 = cleanxppautdat_onepar("src/xppaut/RickerBevertonHolttau3_a.dat")
    sta_data3 = @subset(datatau3, :PointTypeName .== "Stable")
    datatau4 = cleanxppautdat_onepar("src/xppaut/RickerBevertonHolttau4_a.dat")
    sta_data4 = @subset(datatau4, :PointTypeName .== "Stable")
    datatau5 = cleanxppautdat_onepar("src/xppaut/RickerBevertonHolttau5_a.dat")
    sta_data5 = @subset(datatau5, :PointTypeName .== "Stable")
    plot(sta_data0.a, sta_data0.N1, lw=2, linestyle=:solid, color=:blue, label="τ=0")
    plot!(sta_data1.a, sta_data1.N1, lw=2, linestyle=:solid, color=:brown, label="τ=1")
    plot!(sta_data2.a, sta_data2.N1, lw=2, linestyle=:solid, color=:orange, label="τ=2")
    plot!(sta_data3.a, sta_data3.N1, lw=2, linestyle=:solid, color=:purple, label="τ=3")
    plot!(sta_data4.a, sta_data4.N1, lw=2, linestyle=:solid, color=:red, label="τ=4")
    plot!(sta_data5.a, sta_data5.N1, lw=2, linestyle=:solid, color=:green, label="τ=5")
    xlims!(0.0, 500.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("N")
    title!("Ricker Beverton Holt model")
    savefig(joinpath(abpath(), "figs/RickerBevertonHolt_a.pdf"))
end