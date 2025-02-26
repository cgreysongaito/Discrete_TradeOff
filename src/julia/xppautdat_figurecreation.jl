using DelimitedFiles
using DataFrames, DataFramesMeta
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
    data=cleanxppautdat("src/xppaut/testdata.dat")
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
end
