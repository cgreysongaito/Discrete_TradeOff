#Models
function BevertonHolt_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
        return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (a-b*exp(-K*(τ+1)))*(p^(τ+1))*Ndata[t-τ]
end

function BevertonHolt_modelII(Ndata, t, para)
    @unpack α,β,a,b,K,p,D,C,τ = para
        return (Ndata[t]  / (1 + α + β* Ndata[t] )) + (D*((a-b*exp(-K*(τ+1)))*Ndata[t-τ]))/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*(a-b*exp(-K*(τ+1)))*Ndata[t-τ])
end

function RickerConstant_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(τ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + g*(p^(τ+1))*Ndata[t-τ]
end

function RickerConstant_wofec_model(Ndata, t, para; optτ::Int64=0)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(optτ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + g*(p^(τ+1))*Ndata[t-τ]
end

function RickerBeverton_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,D,C,τ = para
    g=a-b*exp(-K*(τ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + (D/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*g*Ndata[t-τ]))*g*Ndata[t-τ]
end

function RickerLeslie_τ0_model(Ndata, t, para)
    @unpack α,β,a,b,K,p,τ = para
    g=a-b*exp(-K*(τ+1))
        return (Ndata[t] * exp(-α-β*Ndata[t])) + g*exp(-α-β*g*Ndata[t])*Ndata[t]
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

function orbitdiagrams(model_func, paraval::String, defaultpar::Union{BevHoltPar, RickerPar}, finalts::Int64; upperval::Union{Float64,Int64}=1.0)
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
        timeseries = model_recursion(10.0, 1000000, local_par, model_func)
        dataN[i] = timeseries[end-finalts:end]
    end
    return [range,dataN]
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
