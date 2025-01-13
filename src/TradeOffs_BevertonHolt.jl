include(packages.jl)

function BevertonHolt_model(para)
        N[t] = N[t-1]  / (1 + α + β* N[t-1] ) + (a-b*exp(-K(τ+1)))*(p^(τ+1))*N[t-τ]
    return N
end

function BevertonHolt_recursion(N0,time,para)
    N = zeros(time)
    N[1] = N0
    for t in 2:time
        N[t] = BevertonHolt_model(N[t-1],para)
    end
    return N
end
#TODO need to figure out history