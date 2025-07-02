#Models
function BevertonHolt_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
        return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

function BevertonHolt_modelII(Ndata, t, para; optτ::Int64=0)
    @unpack α,β,a,b,K,p,D,C,τ = para
    g=calc_g(para)
    return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (D*g*Ndata[t-τ])/((D*(1+D)^(τ+1))+((((1+D)^(τ+1))-1)*C*g*Ndata[t-τ]))
end

function RickerConstant_model(Ndata, t, para; optτ::Int64=0)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(τ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + g*(p^(τ+1))*Ndata[t-τ]
end

function RickerConstant_wofec_model(Ndata, t, para; optτ::Int64=0)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(optτ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + g*(p^(τ+1))*Ndata[t-τ]
end

function RickerLeslie_τ0_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(τ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + g*exp(-α-β*g*Ndata[t])*Ndata[t]
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

function LeslieMatrixOrbitDiagram(τrange, time, finalts, para, init, lesliematrix)
    data = Vector{Vector{Float64}}(undef, length(τrange))
    @threads for τi in eachindex(τrange)
        timeseries = model_Leslierecursion(τrange[τi], time, para, init, lesliematrix)
        data[τi] = first_elements(timeseries)[end-finalts:end]
    end
    return [τrange,data]
end

#TODO switch out second alpha and beta for D and C
#Parameters for the models
@with_kw mutable struct BevHoltPar
    α::Float64 = 0.1 #death rate of mature? check!
    β::Float64 = 0.3 #intraspecific competition rate of mature? check!
    a::Float64 = 10
    b::Float64 = 200
    K::Float64 = 1.0
    p::Float64 = 0.4
    D::Float64 = 0.1 #death rate of immature? check!
    C::Float64 = 0.1 #intraspecific competition rate of immature? check!
    τ::Int64 = 5
end

@with_kw mutable struct RickerPar
    α::Float64 = 0.1
    β::Float64 = 0.3
    a::Float64 = 10.0
    b::Float64 = 200.0
    K::Float64 = 1.0
    p::Float64 = 0.4
    D::Float64 = 0.1
    C::Float64 = 0.1
    τ::Int64 = 5
end

function calc_m(para)
    @unpack a,b,K,p,τ = para
    m = (a-b*exp(-K*(τ+1)))*(p^(τ+1))
    return m
end

function calc_g(para)
    @unpack a,b,K,τ = para
    g = a-b*exp(-K*(τ+1))
    return g
end

function alowerconstraint(para)
    @unpack b, K, τ = para
    return b*exp(-K*(τ+1))
end

function ahigherconstraint(para)
    @unpack a, b, K, τ, p = para
    return (1/(p^(τ+1)))+b*exp(-K*(τ+1))
end

function phigherconstraint(para)
    @unpack a, b, K, τ = para
    return (1/(a-b*exp(-K*(τ+1))))^(1/(τ+1))
end


function model_recursion(N0::Float64, time::Int64, para, model_func; optτ::Int64=0)
    if !isa(time, Int64)
        error("time variable needs to be Int64")
    end
    if !isa(N0, Float64)
        error("N0 variable needs to be Float64")
    end
    @unpack τ = para
    N = fill(N0, τ+1)
    for t in τ+1:τ+time
        newN = model_func(N, t, para; optτ)
        append!(N, newN)
    end
    return N[τ+1:end]
end

function model_τbifurc(τrange, model_func, par::Union{BevHoltPar, RickerPar}, pval)
    data = zeros(length(τrange))
    @threads for i in eachindex(τrange)
        local_par = deepcopy(par)
        local_par.τ = τrange[i]
        local_par.p = pval
        timeseries = model_recursion(0.1, 500, local_par, model_func)
        data[i] = timeseries[end-50]
    end
    return data
end #using end value (#TODO code max 0 or equilibrium point for changing tau)

function model_τorbit(τrange, model_func, par::Union{BevHoltPar, RickerPar}, finalts)
    dataN = Vector{Vector{Float64}}(undef, length(τrange))
    @threads for i in eachindex(τrange)
        local_par = deepcopy(par)
        local_par.τ = τrange[i]
        timeseries = model_recursion(0.1, 500, local_par, model_func)
        dataN[i] = timeseries[end-finalts:end]
    end
    return [τrange,dataN]
end

function orbitdiagrams(model_func, paraval::String, defaultpar::Union{BevHoltPar, RickerPar}, finalts::Int64; optτ::Int64=0, upperval::Union{Float64,Int64}=1.0)
    if paraval == "a" && model_func==RickerConstant_model
        range=round(alowerconstraint(defaultpar), digits=2)+0.1:0.01:round(ahigherconstraint(defaultpar),digits=2)-0.1
    elseif paraval == "a" 
        range=round(alowerconstraint(defaultpar), digits=2)+0.1:1.0:upperval
    elseif paraval == "p"
        upperlimittest=phigherconstraint(defaultpar)
        if upperlimittest>1.0
            upperlimit=1.0
        else
            upperlimit=upperlimittest
        end
        range=0.0:0.001:upperlimit-0.001
    elseif paraval=="α"
        range=0.01:0.01:upperval
    elseif paraval=="β"
        range=0.001:0.001:upperval
    elseif paraval=="τ"
        range=2:1:upperval
    else
        error("paraval should be either a, p, α, β, or τ")
    end
    dataN = Vector{Vector{Float64}}(undef, length(range))
    @threads for i in eachindex(range)
        local_par = deepcopy(defaultpar)
        if paraval == "a"
            local_par.a = range[i]
        elseif paraval == "p"
            local_par.p = range[i]
        elseif paraval == "α"
            local_par.α = range[i]
        elseif paraval == "β"
            local_par.β = range[i]
        else 
            local_par.τ = range[i]
        end
        timeseries = model_recursion(10.0, 1000000, local_par, model_func; optτ)
        dataN[i] = timeseries[end-finalts:end]
    end
    return [range,dataN]
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
#root solve BevertonholtBevertonholt equality
function BHBH_existence_check(para)
    @unpack α,β,D,C,a,b,K,τ = para
    g = calc_g(para)
    return (1-(1/(1+α)))-(g/((1+D)^(τ+1)))
end

function Nequi_BHBH(N, para)
    @unpack α,β,D,C,a,b,K,τ = para
    g = calc_g(para)
    #equilibrium point
    return 1-(1/(1+α+β*N))-(D*g)/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*g*N)
end

function findequil_BHBH(para)
    if BHBH_existence_check(para) > 0
        return 0.0 #no interior equilibrium point exists
    else
    return find_zero(N -> Nequi_BHBH(N, para), 0.0001)
    end
end

function NdataBHBH(CDbetaalpha,τrange,val, defaultpara)
    data = zeros(length(τrange))
    for i in eachindex(τrange)
        local_par = deepcopy(defaultpara)
        local_par.τ = τrange[i]
        if CDbetaalpha == "C"
            local_par.C = val
        elseif CDbetaalpha == "D"
            local_par.D = val
        elseif CDbetaalpha == "β"
            local_par.β = val
        elseif CDbetaalpha == "α"
            local_par.α = val
        else
            error("Invalid parameter type specified. Use 'C', 'D', 'β', or 'α'.")
        end
        data[i] = findequil_BHBH(local_par)
    end
    return data
end

#Accessory functions

# function plot_combination(single_vector, vector_of_vectors, pointcolor)
#     # Loop through the elements and plot the points
#     for i in eachindex(single_vector)
#         for y in eachindex(vector_of_vectors[i])
#             scatter!([single_vector[i]], [vector_of_vectors[i][y]],label=false, color=pointcolor )
#         end
#     end
#     plot!()
# end
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

function flattenorbitdata(orbitdata)
    # Loop through the elements and plot the points
    xaxisdata=[]
    yaxisdata=[]
    single_vector=orbitdata[1]
    vector_of_vectors=orbitdata[2]
    for i in eachindex(single_vector)
        for y in eachindex(vector_of_vectors[i])
            push!(xaxisdata,single_vector[i])
            push!(yaxisdata,vector_of_vectors[i][y])
        end
    end
    return [xaxisdata, yaxisdata]
end

function abpath()
    replace(@__DIR__, "src/julia" => "")
end

# function split_by_min_low2ndpar(df)
#     min_low2ndpar = minimum(df.Low2ndpar)
#     idx = findfirst(==(min_low2ndpar), df.Low2ndpar)
#     a_split = df.a[idx]
#     df_lower = @subset(df, :a .< a_split)
#     df_upper = @subset(df, :a .>= a_split)
#     return df_lower, df_upper, a_split, min_low2ndpar
# end

# df_lower, df_upper, a_split, min_low2ndpar = split_by_min_low2ndpar(hp_datatau1)

