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
    maxmin6=maxminNS(6.0, 32.03:0.01:ahigh6, 1000000, 999000)

    df_maxmin2 = DataFrame(a = maxmin2[1], maximum = maxmin2[2], minimum = maxmin2[3])
    CSV.write("data/maxmin2.csv", df_maxmin2)

    df_maxmin4 = DataFrame(a = maxmin4[1], maximum = maxmin4[2], minimum = maxmin4[3])
    CSV.write("data/maxmin4.csv", df_maxmin4)

    df_maxmin6 = DataFrame(a = maxmin6[1], maximum = maxmin6[2], minimum = maxmin6[3])
    CSV.write("data/maxmin6.csv", df_maxmin6)
end

function findRCvalue(data, aval, stable)
    if stable=="stable"
        return @subset(data, isapprox.(:a, aval; rtol=1e-4)).N1[1]
    elseif stable=="unstable"
        return @subset(data, isapprox.(:a, aval; rtol=1e-4)).N1[1]
    elseif stable=="max"
        return @subset(data, isapprox.(:a, aval; rtol=1e-4)).maximum[1]
    else
        error("Invalid stable type. Use 'stable', 'unstable', or 'max'.")
    end
end

findRCvalue(datatau6, 33.0, "stable")

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
    plot!(maxmindata2.a, maxmindata2.maximum, color="#73D055FF", linestyle=:dashdotdot, lw=2, label="")
    plot!(maxmindata2.a, maxmindata2.minimum, color="#73D055FF", linestyle=:dashdotdot,lw=2, label="")
    plot!(maxmindata4.a, maxmindata4.maximum, color="#1F968BFF", linestyle=:dashdotdot,lw=2, label="")
    plot!(maxmindata4.a, maxmindata4.minimum, color="#1F968BFF", linestyle=:dashdotdot,lw=2, label="")
    plot!(maxmindata6.a, maxmindata6.maximum, color="#404788FF", linestyle=:dashdotdot,lw=2, label="")
    plot!(maxmindata6.a, maxmindata6.minimum, color="#404788FF", linestyle=:dashdotdot,lw=2, label="")
    plot!([-1], [0], linestyle=:dashdotdot, color=:black,lw=2, label="N-S max/min")
    plot!([32, 32], [-0.5, findRCvalue(datatau6, 32.0, "stable")], color=:black, linestyle=:dash, lw=1, label="")
    plot!([33, 33], [-0.5, findRCvalue(maxmindata6, 33.0, "max")], color=:black, linestyle=:dash, lw=1, label="")
    xlabel!("a")
    ylabel!("N*(a)")
    ylims!(-0.5, 20.0)
    xlims!(0.0, 38.0)
    savefig(joinpath(abpath(), "figs/RickerConstanttaueven_a.pdf"))
end

let #Time embedding for tau=6
    endtime = 1000000
    finalts = 999000
    timeseries = model_recursion(0.1, endtime, RickerPar(τ=6, p=0.6, a=33.0, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    df = DataFrame(time = 0:endtime, N = timeseries)
    dfxaxis = @subset(df, :time .>= finalts-6 .&& :time .<= endtime-6)
    dfyaxis = @subset(df, :time .>= finalts)
    dfxaxis.row = 1:nrow(dfxaxis)
    dfyaxis.row = 1:nrow(dfyaxis)
    merged_df = leftjoin(dfxaxis, dfyaxis, on=:row, makeunique=true)
    select!(merged_df, Not([:time, :time_1, :row]))
    timeseriess = model_recursion(0.1, endtime, RickerPar(τ=6, p=0.6, a=32.0, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    dfs = DataFrame(time = 0:endtime, N = timeseriess)
    dfxaxiss = @subset(dfs, :time .>= finalts-6 .&& :time .<= endtime-6)
    dfyaxiss = @subset(dfs, :time .>= finalts)
    dfxaxiss.row = 1:nrow(dfxaxiss)
    dfyaxiss.row = 1:nrow(dfyaxiss)
    merged_dfs = leftjoin(dfxaxiss, dfyaxiss, on=:row, makeunique=true)
    select!(merged_dfs, Not([:time, :time_1, :row]))
    scatter(merged_dfs.N, merged_dfs.N_1, color=:black, label="Fixed Point", ms=8, marker=:star5)
    scatter!(merged_df.N, merged_df.N_1, color=:black, label="N-S", ms=3)
    xlabel!("N(t-τ)")
    ylabel!("N(t)")
    title!("\$\\tau=6\$")
    plot!(size=(400, 300),legend=:bottomleft)
    savefig(joinpath(abpath(), "figs/RickerConstanttau6_timeembedding.pdf"))
end

function branchsplitter(data)
    branches = Vector{DataFrame}()
    current_branch = DataFrame(a=Float64[], N1=Float64[])
    threshold = 1e-2  # adjust as needed

    for i in 1:(nrow(data)-1)
        push!(current_branch, (a=data.a[i], N1=data.N1[i]))
        if abs(data.N1[i+1] - data.N1[i]) > threshold
            push!(branches, deepcopy(current_branch))
            empty!(current_branch)
        end
    end
    # Add the last point and branch
    push!(current_branch, (a=data.a[end], N1=data.N1[end]))
    push!(branches, current_branch)
    # Remove branches with only two rows
    branches = filter(df -> nrow(df) > 2, branches)
    return branches
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
    datatau3u3 = branchsplitter(@subset(datatau3_filtered, :PointTypeName .== "Unstable" .&& :a .> 7.00 .&& :LineNum .== 3.0 .&& :N1 .> 0.00))
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
    plot!(datatau3u3[1].a, datatau3u3[1].N1, color=:black, linestyle=:dash, lw=1.5, label="")

    plot!(datatau3u3[2].a, datatau3u3[2].N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u4upper.a, datatau3u4upper.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u4lower.a, datatau3u4lower.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3u0.a, datatau3u0.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau3ul0.a, datatau3ul0.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    ylabel!("N*(a)")
    xlabel!("a")
    xlims!(tau3lowerbound, tau3upperbound + 0.1)
    ylims!(-0.5, 30.0)
    plot!([11.04, 11.04], [-0.5, findRCvalue(datatau3u3[1], 11.04, "unstable")], color=:black, linestyle=:dash, lw=1, label="")
    plot!([10.92, 10.92], [-0.5, findRCvalue(datatau3s3p4upperupper, 10.92, "stable")], color=:black, linestyle=:dash, lw=1, label="")
    title!("\$\\tau=3\$")
    savefig(joinpath(abpath(), "figs/RickerConstanttau3_a.pdf"))
end

let #Time embedding for tau=3
    endtime = 1000000
    finalts = 999000
    timeseriesp4 = model_recursion(0.1, endtime, RickerPar(τ=3, p=0.6, a=10.92, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
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
    timeseries = model_recursion(0.1, endtime, RickerPar(τ=3, p=0.6, a=11.04, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
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
    
    xlabel!("N(t-τ)")
    ylabel!("N(t)")
    # title!("\$\\tau=3\$")
    savefig(joinpath(abpath(), "figs/RickerConstanttau3_timeembedding.pdf"))
end

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
    datatau5u0 = @subset(datatau5_filtered, :PointTypeName .== "Unstable" .&& :LineNum .== 2.0 .&& :N1 .< 0.001)

    plot(datatau5s1.a, datatau5s1.N1, color=:black, label="Stable")
    plot!(datatau5s2pupper.a, datatau5s2pupper.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau5s2plower.a, datatau5s2plower.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau5s2_0.a, datatau5s2_0.N1, color=:black, linestyle=:solid, label="")
    plot!(datatau5u1.a, datatau5u1.N1, color=:black, linestyle=:dash, lw=1.5, label="Unstable")
    plot!(datatau5u2pupper.a, datatau5u2pupper.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau5u2plower.a, datatau5u2plower.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!(datatau5u0.a, datatau5u0.N1, color=:black, linestyle=:dash, lw=1.5, label="")
    plot!([20.2, 20.2], [-0.5, findRCvalue(datatau5s2pupper, 20.2, "stable")], color=:black, linestyle=:dash, lw=1, label="")
    plot!([20.7, 20.7], [-0.5, findRCvalue(datatau5u2pupper, 20.7, "unstable")], color=:black, linestyle=:dash, lw=1, label="")

    ylabel!("N*(a)")
    xlabel!("a")
    xlims!(0.0, tau5upperbound + 0.1)
    ylims!(-0.5, 30.0)
    title!("\$\\tau=5\$")
    savefig(joinpath(abpath(), "figs/RickerConstanttau5_a.pdf"))
end

let #Time embedding for tau=5
    endtime = 1000000
    finalts = 999000
    timeseriesp2 = model_recursion(0.1, endtime, RickerPar(τ=5, p=0.6, a=20.2, α=0.1, β=0.3, K=1.0, b=200), RickerConstant_model)
    dfp2 = DataFrame(time = 0:endtime, N = timeseriesp2)
    dfp2.periodpoint = mod.(dfp2.time, 2) .+ 1
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
    
    xlabel!("N(t-τ)")
    ylabel!("N(t)")
    # title!("\$\\tau=5\$")
    savefig(joinpath(abpath(), "figs/RickerConstanttau5_timeembedding.pdf"))
end

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
    datatau8 = cleanxppautdat_twopar("src/xppaut/RickerConstanttau8.dat")
    bp_datatau8 = unique(@subset(datatau8, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau8 = unique(@subset(datatau8, :PointTypeName .== "HP"), :Low2ndpar)
    plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau8.a, bp_datatau8.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=8")
    plot!(hp_datatau8.a, hp_datatau8.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="Oscillations")
    scatter!([13.0], [0.57205], color=:black, marker=:star5, ms=12, label="")
    scatter!([13.0], [0.35], color=:black, marker=:square, ms=8, label="")
    xlims!(0.0, 40.0)
    ylims!(0.0, 1.0)
    xlabel!("a")
    ylabel!("\${\\tilde{p}}\$")
    savefig(joinpath(abpath(), "figs/RickerConstant_apbifurcation.pdf"))
end

#Tau create and then kill oscillations
let 
    par = RickerPar(a=13.0, p=0.57205, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; optτ=2, upperval=8))
    timeseries = model_recursion(10.0, 1000000, RickerPar(a=13.0, p=0.57205, τ=1.0, α=0.1, β=0.3, b=200, K=1.0), RickerConstant_wofec_model; optτ=2)
    dataN = timeseries[end-50:end]
    p2 = scatter(RCorbitdata[1], log10.(RCorbitdata[2] .+1), color=:black, ms=8,label="")
    scatter!([1.0], [0.0], color=:black, ms=8, label="")
    xlabel!("τ")
    ylabel!("\$log_{10}(N+1)\$")
    xticks!([0, 2, 4, 6, 8])
    xlims!(0.0, 8.5)
    plot!(size=(350, 300))
    savefig(joinpath(abpath(), "figs/RickerConstant_tauorbita.pdf"))
end

#Tau allows equilibrium to exist and then kills stable point
let
    par = RickerPar(a=13.0, p=0.35, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RCorbitdata = flattenorbitdata(orbitdiagrams(RickerConstant_model, "τ", par, 50; upperval=8))
    timeseries = model_recursion(10.0, 1000000, RickerPar(a=13.0, p=0.35, τ=1.0, α=0.1, β=0.3, b=200, K=1.0), RickerConstant_wofec_model; optτ=2)
    dataN = timeseries[end-50:end]
    p2 = scatter(RCorbitdata[1], RCorbitdata[2], color=:black, ms=8, label="")
    scatter!([1.0], [0.0], color=:black, ms=8,label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 8.5)
    xticks!([0, 2, 4, 6, 8])
    plot!(size=(350, 300))
    savefig(joinpath(abpath(), "figs/RickerConstant_tauorbitb.pdf"))
end

#BevertonHoltBevertonHolt model


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
    datatau9 = cleanxppautdat_twopar("src/xppaut/RickerRickertau9_C.dat")
    bp_datatau9 = unique(@subset(datatau9, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau9 = unique(@subset(datatau9, :PointTypeName .== "HP"), :Low2ndpar)

    p1=plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!(bp_datatau9.a, bp_datatau9.Low2ndpar, lw=2, linestyle=:solid, color="#440154FF", label="τ=9")
    plot!(hp_datatau9.a, hp_datatau9.Low2ndpar, lw=2, linestyle=:dash, color="#440154FF", label="")
    scatter!([26.5], [1.0], color=:black, marker=:star5, ms=12, label="")
    scatter!([26.5], [0.1], color=:black, marker=:square, ms=8, label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="N-S")
    xlims!(0.0, 50.0)
    ylims!(0.0, 5.0)
    xlabel!("a")
    ylabel!("C")
    p2=plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!(bp_datatau9.a, bp_datatau9.Low2ndpar, lw=2, linestyle=:solid, color="#440154FF", label="τ=9")
    plot!(hp_datatau9.a, hp_datatau9.Low2ndpar, lw=2, linestyle=:dash, color="#440154FF", label="")
    scatter!([26.5], [0.1], color=:black, marker=:square, ms=8, label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="N-S")
    xlims!(0.0, 50.0)
    ylims!(0.0, 0.3)
    xlabel!("a")
    ylabel!("C")
    title!("Weak Competition")
    plot(p1, p2, layout=(2, 1), size=(500, 600), titlefont=font(12, "Arial"))
    savefig(joinpath(abpath(), "figs/aCbifurcation_RickerRicker.pdf"))
end

#Tau orbit for a and C comparison (RickerRicker)
let 
    τrange = 1:1:35
    par1 = RickerPar(a=26.5,C=1.0, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RLorbitdata1 = flattenorbitdata(LeslieMatrixOrbitDiagram(τrange, 50000, 50, par1, 0.1, LeslieMatrix))
    par2 = RickerPar(a=26.5,C=0.1, τ=2.0, α=0.1, β=0.3, b=200, K=1.0)
    RLorbitdata2 = flattenorbitdata(LeslieMatrixOrbitDiagram(τrange, 50000, 50, par2, 0.1, LeslieMatrix))
    p1 = scatter(RLorbitdata1[1], RLorbitdata1[2], color=:black, ms=5,label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 35.5)
    p2 = scatter(RLorbitdata2[1], RLorbitdata2[2], color=:black, ms=5,label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 35.5)
    plot(p1, p2, size=(450, 600), layout=(2, 1), legend=:none)
    savefig(joinpath(abpath(), "figs/RickerRicker_aC_tauorbit.pdf"))
end

function split_by_min_low2ndpar(df)
    min_low2ndpar = minimum(df.Low2ndpar)
    idx = findfirst(==(min_low2ndpar), df.Low2ndpar)
    a_split = df.a[idx]
    df_lower = @subset(df, :a .< a_split)
    df_upper = @subset(df, :a .>= a_split)
    return df_lower, df_upper, a_split, min_low2ndpar
end

df_lower, df_upper, a_split, min_low2ndpar = split_by_min_low2ndpar(hp_datatau1)


#a and β two par bifurcation figure for RickerRicker
let
    datatau1 = cleanxppautdat_twopar("src/xppaut/RickerRickertau1_beta_b.dat")
    bp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau1 = unique(@subset(datatau1, :PointTypeName .== "HP" .&& :a .<= 70), :Low2ndpar)
    datatau2 = cleanxppautdat_twopar("src/xppaut/RickerRickertau2_beta_b.dat")
    bp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau2 = unique(@subset(datatau2, :PointTypeName .== "HP"), :Low2ndpar)
    datatau3 = cleanxppautdat_twopar("src/xppaut/RickerRickertau3_beta_b.dat")
    bp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau3 = unique(@subset(datatau3, :PointTypeName .== "HP"), :Low2ndpar)
    datatau4 = cleanxppautdat_twopar("src/xppaut/RickerRickertau4_beta_b.dat")
    bp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau4 = unique(@subset(datatau4, :PointTypeName .== "HP"), :Low2ndpar)
    datatau5 = cleanxppautdat_twopar("src/xppaut/RickerRickertau5_beta_b.dat")
    bp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau5 = unique(@subset(datatau5, :PointTypeName .== "HP"), :Low2ndpar)
    datatau9 = cleanxppautdat_twopar("src/xppaut/RickerRickertau9_beta_b.dat")
    bp_datatau9 = unique(@subset(datatau9, :PointTypeName .== "BP"), :Low2ndpar)
    hp_datatau9 = unique(@subset(datatau9, :PointTypeName .== "HP"), :Low2ndpar)
    p1=plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!(bp_datatau9.a, bp_datatau9.Low2ndpar, lw=2, linestyle=:solid, color="#440154FF", label="τ=9")
    plot!(hp_datatau9.a, hp_datatau9.Low2ndpar, lw=2, linestyle=:dash, color="#440154FF", label="")
    scatter!([26.5], [1.0], color=:black, marker=:star5, ms=10, label="")
    scatter!([26.5], [0.1], color=:black, marker=:square, ms=6, label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="N-S")
    xlims!(0.0, 70.0)
    ylims!(0.0, 5.0)
    xlabel!("a")
    ylabel!("β")
    p2=plot(bp_datatau1.a, bp_datatau1.Low2ndpar, lw=2, linestyle=:solid, color="#FDE725FF", label="τ=1")
    plot!(hp_datatau1.a, hp_datatau1.Low2ndpar, lw=2, linestyle=:dash, color="#FDE725FF", label="")
    plot!(bp_datatau2.a, bp_datatau2.Low2ndpar, lw=2, linestyle=:solid, color="#73D055FF", label="τ=2")
    plot!(hp_datatau2.a, hp_datatau2.Low2ndpar, lw=2, linestyle=:dash, color="#73D055FF", label="")
    plot!(bp_datatau3.a, bp_datatau3.Low2ndpar, lw=2, linestyle=:solid, color="#20A387FF", label="τ=3")
    plot!(hp_datatau3.a, hp_datatau3.Low2ndpar, lw=2, linestyle=:dash, color="#20A387FF", label="")
    plot!(bp_datatau4.a, bp_datatau4.Low2ndpar, lw=2, linestyle=:solid, color="#2D708EFF", label="τ=4")
    plot!(hp_datatau4.a, hp_datatau4.Low2ndpar, lw=2, linestyle=:dash, color="#2D708EFF", label="")
    plot!(bp_datatau5.a, bp_datatau5.Low2ndpar, lw=2, linestyle=:solid, color="#453781FF", label="τ=5")
    plot!(hp_datatau5.a, hp_datatau5.Low2ndpar, lw=2, linestyle=:dash, color="#453781FF", label="")
    plot!(bp_datatau9.a, bp_datatau9.Low2ndpar, lw=2, linestyle=:solid, color="#440154FF", label="τ=9")
    plot!(hp_datatau9.a, hp_datatau9.Low2ndpar, lw=2, linestyle=:dash, color="#440154FF", label="")
    scatter!([26.5], [0.1], color=:black, marker=:square, ms=6, label="")
    plot!([NaN], [NaN], lw=2, linestyle=:solid, color=:black, label="Transcritical")
    plot!([NaN], [NaN], lw=2, linestyle=:dash, color=:black, label="N-S")
    xlims!(0.0, 70.0)
    ylims!(0.0, 0.3)
    xlabel!("a")
    ylabel!("β")
    title!("Weak Competition")
    plot!(legend=:topright)
    plot(p1, p2, layout=(2, 1), size=(500, 600), titlefont=font(12, "Arial"))
    savefig(joinpath(abpath(), "figs/abetabifurcation_RickerRicker.pdf"))
end

#Tau orbit for a and β comparison (RickerRicker)
let 
    τrange = 1:1:35
    par1 = RickerPar(a=26.5, β=1.0, C=0.3, τ=2.0, α=0.1, b=200, K=1.0)
    RLorbitdata1 = flattenorbitdata(LeslieMatrixOrbitDiagram(τrange, 50000, 50, par1, 0.1, LeslieMatrix))
    par2 = RickerPar(a=26.5, β=0.1, C=0.3, τ=2.0, α=0.1, b=200, K=1.0)
    RLorbitdata2 = flattenorbitdata(LeslieMatrixOrbitDiagram(τrange, 50000, 50, par2, 0.1, LeslieMatrix))
    p1 = scatter(RLorbitdata1[1], RLorbitdata1[2], color=:black, ms=5,label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 35.5)
    p2 = scatter(RLorbitdata2[1], RLorbitdata2[2], color=:black, ms=5,label="")
    xlabel!("τ")
    ylabel!("N")
    xlims!(0.0, 35.5)
    plot(p1, p2, layout=(2, 1), size=(450, 600), titlefont=font(12, "Arial"))
    savefig(joinpath(abpath(), "figs/RickerRicker_ab_tauorbit.pdf"))
end

#Appendix

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