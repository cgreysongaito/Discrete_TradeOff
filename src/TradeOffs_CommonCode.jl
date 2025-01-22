#Parameters for the models
@with_kw mutable struct BevHoltPar
    α::Float64 = 0.1 #death rate of mature? check!
    β::Float64 = 0.3 #intraspecific competition rate of mature? check!
    a::Float64 = 10
    b::Float64 = 200
    K::Float64 = 1
    p::Float64 = 0.4
    D::Float64 = 0.1 #death rate of immature? check!
    C::Float64 = 0.1 #intraspecific competition rate of immature? check!
    τ::Int = 5
end

@with_kw mutable struct RickerPar
    r::Float64 = 0.1
    β::Float64 = 0.3
    a::Float64 = 10
    b::Float64 = 200
    K::Float64 = 1
    p::Float64 = 0.4
    D::Float64 = 0.1
    C::Float64 = 0.1
    τ::Int = 5
end


function model_recursion(N0, time, para, model_func)
    @unpack τ = para
    N = fill(N0, τ+1)
    for t in τ+1:τ+time
        newN = model_func(N, t, para)
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
end
# function model_τbifurc(τrange, model_func, par)
#     data= zeros(length(τrange))
#     @threads for i in eachindex(τrange)
#         timeseries = model_recursion(0.1,500, BevHoltPar(τ=τrange[i], p=pval), model_func)
#         data[i] = timeseries[end-50]
#     end
#     return data
# end

