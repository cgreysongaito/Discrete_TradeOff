include("packages.jl")
# default(grid=false, linewidth=3, tickfontsize=12, legendfontsize=10, guidefontsize=15)
include("TradeOffs_CommonCode.jl")
using CSV
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
#Calculate max and min for Neimark-Sacker
function maxminNS(τval, arange, time, trans)
    maxdata=zeros(length(arange))
    mindata=zeros(length(arange))
    @threads for ai in eachindex(arange)
        timeseries = model_recursion(0.1, time, RickerPar(τ=τval, a=arange[ai], p=0.6, α=0.1, β=0.3, b=200, K=1.0), RickerConstant_model)[trans:end]
        maxdata[ai] = maximum(timeseries)
        mindata[ai] = minimum(timeseries)
    end
    return arange, maxdata, mindata
end

let #create max min N-S data for even tau Ricker constant
    ahigh2=ahigherconstraint(RickerPar(p=0.6, τ=2.0, α=0.1, β=0.3, b=200, K=1.0))
    maxmin2=maxminNS(2.0, 14.39:0.01:ahigh2, 1000000, 999000)
    ahigh4=ahigherconstraint(RickerPar(p=0.6, τ=4.0, α=0.1, β=0.3, b=200, K=1.0))
    maxmin4=maxminNS(4.0, 12.99:0.01:ahigh4, 1000000, 999000)
    ahigh6=ahigherconstraint(RickerPar(p=0.6, τ=6.0, α=0.1, β=0.3, b=200, K=1.0))
    maxmin6=maxminNS(6.0, 32.05:0.01:ahigh6, 1000000, 999000)

    df_maxmin2 = DataFrame(arange = maxmin2[1], maximum = maxmin2[2], minimum = maxmin2[3])
    CSV.write("data/maxmin2.csv", df_maxmin2)

    df_maxmin4 = DataFrame(arange = maxmin4[1], maximum = maxmin4[2], minimum = maxmin4[3])
    CSV.write("data/maxmin4.csv", df_maxmin4)

    df_maxmin6 = DataFrame(arange = maxmin6[1], maximum = maxmin6[2], minimum = maxmin6[3])
    CSV.write("data/maxmin6.csv", df_maxmin6)
end

let 
    timeseries = model_recursion(0.1, 1000000, RickerPar(τ=6.0, a=35.0, p=0.6, α=0.1, β=0.3, b=200, K=1.0), RickerConstant_model)[990000:end]
    println("Max: ", maximum(timeseries), " Min: ", minimum(timeseries))
    plot(timeseries, color=:black, label="τ=6.0, a=32.05", linewidth=2)
    xlabel!("Time")
    ylabel!("N")
end


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
    ahigh6=ahigherconstraint(RickerPar(p=0.6, τ=6.0, α=0.1, β=0.3, b=200, K=1.0))
    maxmindata2=CSV.read("data/maxmin2.csv", DataFrame)
    maxmindata4=CSV.read("data/maxmin4.csv", DataFrame)
    maxmindata6=CSV.read("data/maxmin6.csv", DataFrame)

    plot(datatau2s.a, datatau2s.N1, color="#73D055FF", label="τ=2.0")
    plot!(datatau2ua.a, datatau2ua.N1, color="#73D055FF", linestyle=:dash, label="")
    plot!(datatau4s.a, datatau4s.N1, color="#1F968BFF", label="τ=4.0")
    plot!(datatau4ua.a, datatau4ua.N1, color="#1F968BFF", linestyle=:dash, label="")
    plot!(datatau6s.a, datatau6s.N1, color="#404788FF", label="τ=6.0")
    plot!(datatau6ua.a, datatau6ua.N1, color="#404788FF", linestyle=:dash, label="")
    plot!([-1], [0], linestyle=:solid, color=:black, label="Stable")
    plot!([-1], [0], linestyle=:dash, color=:black, label="Unstable")
    plot!(maxmindata2.arange, maxmindata2.maximum, color=:black, linestyle=:dot, label="Max")
    plot!(maxmindata2.arange, maxmindata2.minimum, color=:black, linestyle=:dot, label="Min")
    plot!(maxmindata4.arange, maxmindata4.maximum, color=:black, linestyle=:dot, label="Max")
    plot!(maxmindata4.arange, maxmindata4.minimum, color=:black, linestyle=:dot, label="Min")
    
    plot!(maxmindata6.arange, maxmindata6.maximum, color=:black, linestyle=:dot, label="Max")
    plot!(maxmindata6.arange, maxmindata6.minimum, color=:black, linestyle=:dot, label="Min")
    xlabel!("a")
    ylabel!("N*")
    ylims!(-0.5, 50.0)
    xlims!(0.0, 38.0)
    # println(maximum(datatau4s.a))
    # savefig(joinpath(abpath(), "figs/RickerConstanttaueven_a.pdf"))
end

let #Time embedding for tau=6
    endtime = 1000000
    finalts = 999000
    timeseries = model_recursion(0.1, endtime, RickerPar(τ=6, p=0.6, a=33, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    df = DataFrame(time = 0:endtime, N = timeseries)
    dfxaxis = @subset(df, :time .>= finalts-6 .&& :time .<= endtime-6)
    dfyaxis = @subset(df, :time .>= finalts)
    dfxaxis.row = 1:nrow(dfxaxis)
    dfyaxis.row = 1:nrow(dfyaxis)
    merged_df = leftjoin(dfxaxis, dfyaxis, on=:row, makeunique=true)
    select!(merged_df, Not([:time, :time_1, :row]))
    timeseriess = model_recursion(0.1, endtime, RickerPar(τ=6, p=0.6, a=32, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    dfs = DataFrame(time = 0:endtime, N = timeseriess)
    dfxaxiss = @subset(dfs, :time .>= finalts-6 .&& :time .<= endtime-6)
    dfyaxiss = @subset(dfs, :time .>= finalts)
    dfxaxiss.row = 1:nrow(dfxaxiss)
    dfyaxiss.row = 1:nrow(dfyaxiss)
    merged_dfs = leftjoin(dfxaxiss, dfyaxiss, on=:row, makeunique=true)
    select!(merged_dfs, Not([:time, :time_1, :row]))

    scatter(merged_df.N, merged_df.N_1, color=:black, label="a=33", ms=3)
    scatter!(merged_dfs.N, merged_dfs.N_1, color=:black, label="a=32", ms=8, marker=:star5)
    xlabel!("N(t-6)")
    ylabel!("N(t)")

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
    # println(maximum(datatau3s3p4lowerupper.a))
end

let #Time embedding for tau=3
    endtime = 1000000
    finalts = 999000
    timeseriesp4 = model_recursion(0.1, endtime, RickerPar(τ=3, p=0.6, a=10.96, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    dfp4 = DataFrame(time = 0:endtime, N = timeseriesp4)
    dfp4.periodpoint = mod.(dfp4.time, 4) .+ 1
    dfp4.periodpoint = map(x -> x == 1 ? 3 : x == 2 ? 4 : x == 3 ? 1 : x == 4 ? 2 : x, dfp4.periodpoint)
    dfxaxisp4 = @subset(dfp4, :time .>= finalts-3 .&& :time .<= endtime-3)
    select!(dfxaxisp4, Not(:periodpoint))
    dfyaxisp4 = @subset(dfp4, :time .>= finalts)
    dfxaxisp4.row = 1:nrow(dfxaxisp4)
    dfyaxisp4.row = 1:nrow(dfyaxisp4)
    merged_dfp4 = leftjoin(dfxaxisp4, dfyaxisp4, on=:row, makeunique=true)
    select!(merged_dfp4, Not([:time, :time_1, :row]))
    df1p4 = @subset(merged_dfp4, :periodpoint .== 1)
    df2p4 = @subset(merged_dfp4, :periodpoint .== 2)
    df3p4 = @subset(merged_dfp4, :periodpoint .== 3)
    df4p4 = @subset(merged_dfp4, :periodpoint .== 4)
    timeseries = model_recursion(0.1, endtime, RickerPar(τ=3, p=0.6, a=10.99, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    df = DataFrame(time = 0:endtime, N = timeseries)
    df.periodpoint = mod.(df.time, 4) .+ 1
    dfxaxis = @subset(df, :time .>= finalts-3 .&& :time .<= endtime-3)
    select!(dfxaxis, Not(:periodpoint))
    dfyaxis = @subset(df, :time .>= finalts)
    dfxaxis.row = 1:nrow(dfxaxis)
    dfyaxis.row = 1:nrow(dfyaxis)
    merged_df = leftjoin(dfxaxis, dfyaxis, on=:row, makeunique=true)
    select!(merged_df, Not([:time, :time_1, :row]))
    df1 = @subset(merged_df, :periodpoint .== 1)
    df2 = @subset(merged_df, :periodpoint .== 2)
    df3 = @subset(merged_df, :periodpoint .== 3)
    df4 = @subset(merged_df, :periodpoint .== 4)
    scatter(df1.N, df1.N_1, color="#FDE725FF", label="1", ms=3)
    scatter!(df2.N, df2.N_1, color="#73D055FF", label="2", ms=3)
    scatter!(df3.N, df3.N_1, color="#238A8DFF", label="3", ms=3)
    scatter!(df4.N, df4.N_1, color="#440154FF", label="4", ms=3)

    scatter!(df1p4.N, df1p4.N_1, color="#FDE725FF", label="", ms=8, marker=:star5)
    scatter!(df2p4.N, df2p4.N_1, color="#73D055FF", label="", ms=8, marker=:star5)
    scatter!(df3p4.N, df3p4.N_1, color="#238A8DFF", label="", ms=8, marker=:star5)
    scatter!(df4p4.N, df4p4.N_1, color="#440154FF", label="", ms=8, marker=:star5)
    scatter!([NaN], [NaN], colour=:black,marker=:star5, label="Period Four")
    scatter!([NaN], [NaN], colour=:black, marker=:circle, label="N-S")
    
    xlabel!("N(t-3)")
    ylabel!("N(t)")

end

# let #Time embedding for tau=3 when period N-S
#     endtime = 1000000
#     finalts = 999000
#     timeseries = model_recursion(0.1, endtime, RickerPar(τ=3, p=0.6, a=10.99, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
#     df = DataFrame(time = 0:endtime, N = timeseries)
#     df.periodpoint = mod.(df.time, 4) .+ 1
#     dfxaxis = @subset(df, :time .>= finalts-3 .&& :time .<= endtime-3)
#     select!(dfxaxis, Not(:periodpoint))
#     dfyaxis = @subset(df, :time .>= finalts)
#     dfxaxis.row = 1:nrow(dfxaxis)
#     dfyaxis.row = 1:nrow(dfyaxis)
#     merged_df = leftjoin(dfxaxis, dfyaxis, on=:row, makeunique=true)
#     select!(merged_df, Not([:time, :time_1, :row]))
#     df1 = @subset(merged_df, :periodpoint .== 1)
#     df2 = @subset(merged_df, :periodpoint .== 2)
#     df3 = @subset(merged_df, :periodpoint .== 3)
#     df4 = @subset(merged_df, :periodpoint .== 4)
#     scatter(df1.N, df1.N_1, color="#FDE725FF", label="1", ms=3)
#     scatter!(df2.N, df2.N_1, color="#73D055FF", label="2", ms=3)
#     scatter!(df3.N, df3.N_1, color="#238A8DFF", label="3", ms=3)
#     scatter!(df4.N, df4.N_1, color="#440154FF", label="4", ms=3)
#     xlabel!("N(t-3)")
#     ylabel!("N(t)")
#     # savefig(joinpath(abpath(), "figs/timeemedding_a108_Rickerconstant.pdf"))
# end

let #tau=5 (odd)
    datatau5 = cleanxppautdat_onepar("src/xppaut/RickerConstanttau5_a.dat")
    tau5lowerbound = alowerconstraint(RickerPar(p=0.6, τ=5.0, α=0.1, β=0.3, b=200, K=1.0))
    tau5upperbound = ahigherconstraint(RickerPar(p=0.6, τ=5.0, α=0.1, β=0.3, b=200, K=1.0))
    datatau5_filtered = sort(@subset(datatau5, :a .> tau5lowerbound), :a)
    datatau5s1 = @subset(datatau5_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 1.0)
    datatau5s2pupper = unique(@subset(datatau5_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 2.0 .&& :N1 .> maximum(datatau5s1.N1)), :a)
    datatau5s2plower = unique(@subset(datatau5_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 2.0 .&& :N1 .< maximum(datatau5s1.N1) .&& :N1 .> 0.1), :a)
    datatau5s2_0 = @subset(datatau5_filtered, :PointTypeName .== "Stable" .&& :LineNum .== 2.0 .&& :N1 .< 1.00)

    datatau5u1 = unique(@subset(datatau5_filtered, :PointTypeName .== "Unstable" .&& :a .> 15 .&& :LineNum .== 1.0 .&& :N1 .> 0.00), :a)
    datatau5u2pupper = @subset(datatau5_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .> maximum(datatau5s1.N1) .&& :a .> maximum(datatau5s1.a) + 0.1)
    datatau5u2plower = @subset(datatau5_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .< maximum(datatau5s1.N1) .&& :N1 .> 0 .&& :a .> maximum(datatau5s1.a) + 0.1)
    # datatau3u3 = unique(@subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 3.0 .&& :N1 .> 0.00), :a)
    # datatau3u4upper = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 4.0 .&& :N1 .> maximum(datatau3s3p4lowerlower.N1))
    # datatau3u4lower = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 4.0 .&& :N1 .< maximum(datatau3s3p4lowerlower.N1) .&& :N1 .> 0.001)
    datatau5u0 = @subset(datatau5_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .< 0.001)
    # datatau3ul0 = @subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 1.0 .&& :N1 .< 0.001)

    plot(datatau5s1.a, datatau5s1.N1, color=:black, label="Stable")
    plot!(datatau5s2pupper.a, datatau5s2pupper.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau5s2plower.a, datatau5s2plower.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau5s2_0.a, datatau5s2_0.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau5u1.a, datatau5u1.N1, color=:black, linestyle=:dash, lw=1.5, label="Unstable")
    plot!(datatau5u2pupper.a, datatau5u2pupper.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau5u2plower.a, datatau5u2plower.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau5u0.a, datatau5u0.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    ylabel!("N*")
    xlabel!("a")
    xlims!(0.0, tau5upperbound + 0.1)
    ylims!(-0.5, 30.0)
    # println(maximum(datatau5s2plower.a))
end

let #Time embedding for tau=5
    endtime = 1000000
    finalts = 999000
    timeseriesp2 = model_recursion(0.1, endtime, RickerPar(τ=5, p=0.6, a=20.5, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    dfp2 = DataFrame(time = 0:endtime, N = timeseriesp2)
    dfp2.periodpoint = mod.(dfp2.time, 2) .+ 1
    # dfp2.periodpoint = map(x -> x == 1 ? 3 : x == 2 ? 4 : x == 3 ? 1 : x == 4 ? 2 : x, dfp4.periodpoint)
    dfxaxisp2 = @subset(dfp2, :time .>= finalts-5 .&& :time .<= endtime-5)
    select!(dfxaxisp2, Not(:periodpoint))
    dfyaxisp2 = @subset(dfp2, :time .>= finalts)
    dfxaxisp2.row = 1:nrow(dfxaxisp2)
    dfyaxisp2.row = 1:nrow(dfyaxisp2)
    merged_dfp2 = leftjoin(dfxaxisp2, dfyaxisp2, on=:row, makeunique=true)
    select!(merged_dfp2, Not([:time, :time_1, :row]))
    df1p2 = @subset(merged_dfp2, :periodpoint .== 1)
    df2p2 = @subset(merged_dfp2, :periodpoint .== 2)
    timeseries = model_recursion(0.1, endtime, RickerPar(τ=5, p=0.6, a=20.7, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    df = DataFrame(time = 0:endtime, N = timeseries)
    df.periodpoint = mod.(df.time, 2) .+ 1
    dfxaxis = @subset(df, :time .>= finalts-5 .&& :time .<= endtime-5)
    select!(dfxaxis, Not(:periodpoint))
    dfyaxis = @subset(df, :time .>= finalts)
    dfxaxis.row = 1:nrow(dfxaxis)
    dfyaxis.row = 1:nrow(dfyaxis)
    merged_df = leftjoin(dfxaxis, dfyaxis, on=:row, makeunique=true)
    select!(merged_df, Not([:time, :time_1, :row]))
    df1 = @subset(merged_df, :periodpoint .== 1)
    df2 = @subset(merged_df, :periodpoint .== 2)
    scatter(df1.N, df1.N_1, color="#FDE725FF", label="1", ms=3)
    scatter!(df2.N, df2.N_1, color="#73D055FF", label="2", ms=3)

    scatter!(df1p2.N, df1p2.N_1, color="#FDE725FF", label="", ms=8, marker=:star5)
    scatter!(df2p2.N, df2p2.N_1, color="#73D055FF", label="", ms=8, marker=:star5)
    scatter!([NaN], [NaN], colour=:black,marker=:star5, label="Period Two")
    scatter!([NaN], [NaN], colour=:black, marker=:circle, label="N-S")
    
    xlabel!("N(t-5)")
    ylabel!("N(t)")

end

# let
#     alow = alowerconstraint(RickerPar(p=0.3, τ=2.0, α=0.1, β=0.3, b=200, K=1.0))
#     arange = alow:0.1:15.0
#     fillbottom = zeros(length(arange))
#     data = cleanxppautdat_twopar("src/xppaut/RickerConstanttau2.dat")
#     bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
#     hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
#     phigh = [phigherconstraint(RickerPar(a=aval, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
#     plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
#     plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
#     plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
#     xlabel!("a")
#     ylabel!("p")
#     xlims!(9.8, 14.5)
#     ylims!(0.0, 1.0)
#     title!("τ = 2.0")
# end

# let
#     alow = alowerconstraint(RickerPar(p=0.3, τ=3.0, α=0.1, β=0.3, b=200, K=1.0))
#     arange = alow:0.1:15.0
#     fillbottom = zeros(length(arange))
#     data = cleanxppautdat("src/xppaut/RickerConstanttau3woPD.dat")
#     bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
#     hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
#     phigh = [phigherconstraint(RickerPar(a=aval, τ=3.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
#     pddata = perioddoublecurve(4.3:0.1:15.0, 3)
#     # return pddata
#     plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
#     plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
#     plot!(pddata[:, 1], pddata[:, 2], label="Period doubling", lw=2)
#     plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
#     xlabel!("a")
#     ylabel!("p")
#     xlims!(0.0, 15.0)
#     ylims!(0.0, 1.0)
#     title!("τ = 3.0")
#     savefig(joinpath(abpath(), "figs/RickerConstanttau3wPD.pdf"))
# end

# let
#     alow = alowerconstraint(RickerPar(p=0.3, τ=4.0, α=0.1, β=0.3, b=200, K=1.0))
#     arange = alow:0.1:15.0
#     fillbottom = zeros(length(arange))
#     data = cleanxppautdat("src/xppaut/RickerConstanttau4.dat")
#     bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
#     hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
#     phigh = [phigherconstraint(RickerPar(a=aval, τ=4.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
#     plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
#     plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
#     plot!(pddata[:, 1], pddata[:, 2], label="Period doubling", lw=2)
#     plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
#     xlabel!("a")
#     ylabel!("p")
#     xlims!(0.0, 15.0)
#     ylims!(0.0, 1.0)
#     title!("τ = 4.0")
# end

# let
#     alow = alowerconstraint(RickerPar(p=0.3, τ=5.0, α=0.1, β=0.3, b=200, K=1.0))
#     arange = alow:0.1:22.0
#     fillbottom = zeros(length(arange))
#     data = cleanxppautdat("src/xppaut/RickerConstanttau5.dat")
#     bp_data = unique(@subset(data, :PointTypeName .== "BP"), :Lowp)
#     hp_data = unique(@subset(data, :PointTypeName .== "HP"), :Lowp)
#     pddata = perioddoublecurve(5.6:0.1:22.0, 5)
#     phigh = [phigherconstraint(RickerPar(a=aval, τ=5.0, α=0.1, β=0.3, b=200, K=1.0)) for aval in arange]
#     plot(bp_data.a, bp_data.Lowp, label="Transcritical", lw=2)
#     plot!(hp_data.a, hp_data.Lowp, label="Neimarck-Sacker", lw=2)
#     plot!(pddata[:, 1], pddata[:, 2], label="Period doubling", lw=2)
#     plot!(arange, fillbottom, fillrange=phigh, fillalpha=0.2, c=1, label="Parameter space")
#     xlabel!("a")
#     ylabel!("p")
#     xlims!(0.0, 22.0)
#     ylims!(0.0, 1.0)
#     title!("τ = 5.0")
# end

let #dimension 2 (a & p) bifurcation diagram of Ricker Constant 
    datatau1 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau1.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau2.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau3.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau4.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    # datatau5 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau5.dat")
    # bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    # hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    # datatau6 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau6.dat")
    # bp_datatau6 = unique(@subset(datatau6, :PointTypeName .== "BP"), :Low2ndpar)
    # hp_datatau6 = unique(@subset(datatau6, :PointTypeName .== "HP"), :Low2ndpar)
    datatau8 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau8.dat")
    bp_datatau8 = unique(@subset(datatau8, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau8 = unique(@subset(datatau8, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color=:black, label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color=:black, label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color=:blue, label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color=:blue, label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color=:purple, label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color=:purple, label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color=:red, label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color=:red, label="")
    # plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color=:green, label="τ=5")
    # plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color=:green, label="")
    # plot!(bp_datatau6.a, bp_datatau6.Low2ndpar, lw=2, linestyle=:solid, color=:orange, label="τ=6")
    # plot!(hp_datatau6.a, hp_datatau6.Low2ndpar, lw=2, linestyle=:dash, color=:orange, label="")
    plot!(bp_datatau8.a, bp_datatau8.Low2ndpar, lw=2, linestyle=:solid, color=:brown, label="τ=8")
    plot!(hp_datatau8.a, hp_datatau8.Low2ndpar, lw=2, linestyle=:dash, color=:brown, label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="Oscillations")
    scatter!([13.0], [0.57205], color=:black, marker=:star5, ms=12, label="")
    scatter!([13.0], [0.35], color=:black, marker=:square, ms=8, label="")
    xlims!(0.0, 40.0)
    ylims!(0.0, 1.0)
    xlabel!("a")
    ylabel!("p")
    plot!(legendfontsize=10)
    savefig(joinpath(abpath(), "figs/RickerConstant_apbifurcation.pdf"))
end

#Tau create and then kill oscillations
let 
    par = RickerPar(a=13.0, p=0.57205, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; optτ=2, upperval=8))
    RCorbitdata_wofec = flattenorbitdata(orbitdiagrams(RickerConstant_wofec_model, "τ", par, 50; optτ=2, upperval=8))
    timeseries = model_recursion(10.0, 1000000, RickerPar(a=13.0, p=0.57205, τ=1.0, α=0.1, β=0.3, b=200, K=1.0), RickerConstant_wofec_model; optτ=2)
    dataN = timeseries[end-50:end]
    p2 = scatter(RCorbitdata[1], log10.(RCorbitdata[2] .+1), color=:black, ms=8,label="With fecundity\nbenefit")
    scatter!(RCorbitdata_wofec[1], log10.(RCorbitdata_wofec[2] .+1), color=:blue, ms=5, markerstrokecolor=:blue, label="Without fecundity\nbenefit")
    scatter!([1.0], [0.0], color=:black, ms=8, label="")
    scatter!(repeat([1],51), log10.(dataN .+1), color=:blue, ms=5, markerstrokecolor=:blue, label="")
    xlabel!("τ")
    ylabel!("log10(N+1)")
    xticks!([0, 2, 4, 6, 8])
    xlims!(0.0, 8.5)
    plot!(size=(350, 300))
    savefig(joinpath(abpath(), "figs/RickerConstant_tauorbita.pdf"))
end

#Tau allows equilibrium to exist and then kills stable point
let
    par = RickerPar(a=13.0, p=0.35, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; upperval=8))
    RCorbitdata_wofec = flattenorbitdata(orbitdiagrams(RickerConstant_wofec_model, "τ", par, 50; optτ=2, upperval=8))
    timeseries = model_recursion(10.0, 1000000, RickerPar(a=13.0, p=0.35, τ=1.0, α=0.1, β=0.3, b=200, K=1.0), RickerConstant_wofec_model; optτ=2)
    dataN = timeseries[end-50:end]
    p2 = scatter(RCorbitdata[1], RCorbitdata[2], color=:black, ms=8, label="With fecundity\nbenefit")
    scatter!(RCorbitdata_wofec[1], log10.(RCorbitdata_wofec[2] .+1), color=:blue, ms=5, markerstrokecolor=:blue, label="Without fecundity\nbenefit")
    scatter!([1.0], [0.0], color=:black, ms=8,label="")
    scatter!(repeat([1],51), log10.(dataN .+1), color=:blue, ms=5, markerstrokecolor=:blue, label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 8.5)
    xticks!([0, 2, 4, 6, 8])
    plot!(size=(350, 300))
    savefig(joinpath(abpath(), "figs/RickerConstant_tauorbitb.pdf"))
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
    # savefig(joinpath(abpath(), "figs/BevertonHoltBevertonHolt_a.pdf"))
end

#RickerRicker
#a and C two par bifurcation figure

let
    datatau1 = cleanxppautdat_twopar("src/xppaut/RickerRickertau1_C.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat_twopar("src/xppaut/RickerRickertau2_C.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat_twopar("src/xppaut/RickerRickertau3_C.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat_twopar("src/xppaut/RickerRickertau4_C.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat_twopar("src/xppaut/RickerRickertau5_C.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="N-S")
    xlims!(0.0, 50.0)
    ylims!(0.0, 10.0)
    xlabel!("a")
    ylabel!("C")
    savefig(joinpath(abpath(), "figs/aCbifurcation_RickerRicker.pdf"))
end

let
    datatau1 = cleanxppautdat_twopar("src/xppaut/RickerRickertau1_C.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP"), :Low2ndpar)
    datatau2 = cleanxppautdat_twopar("src/xppaut/RickerRickertau2_C.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat_twopar("src/xppaut/RickerRickertau3_C.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat_twopar("src/xppaut/RickerRickertau4_C.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat_twopar("src/xppaut/RickerRickertau5_C.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="N-S")
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