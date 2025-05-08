include("packages.jl")
include("TradeOffs_CommonCode.jl")

function BHMature(N, par)
    @unpack α, β = par
    return 1/(1+α+β*N)
end

function BHImmature(N, par)
    @unpack D, C, τ = par
    g=calc_g(par)
    return (D/((D*(1+D)^(τ+1))+(((1+D)^(τ+1))-1)*C*g*N))
end

let 
    Nrange = 0:0.1:10
    BHMdata = [BHMature(N, BevHoltPar()) for N in Nrange]
    plot(Nrange, BHMdata)
    xlabel!("N")
    ylabel!("p(Nt)")
end

let 
    Nrange = 0:0.1:10
    BHImdatatau0 = [BHImmature(N, BevHoltPar(τ=0, a=75.0)) for N in Nrange]
    BHImdatatau1 = [BHImmature(N, BevHoltPar(τ=1, a=75.0)) for N in Nrange]
    BHImdatatau2 = [BHImmature(N, BevHoltPar(τ=2, a=75.0)) for N in Nrange]
    BHImdatatau3 = [BHImmature(N, BevHoltPar(τ=3, a=75.0)) for N in Nrange]
    plot(Nrange, BHImdatatau0, label="τ=0")
    plot!(Nrange, BHImdatatau1, label="τ=1")
    plot!(Nrange, BHImdatatau2, label="τ=2")
    plot!(Nrange, BHImdatatau3, label="τ=3")
    xlabel!("N")
    ylabel!("p(Nt-τ)")
end

calc_g(BevHoltPar(τ=0, a=75.0))

function RMature(N, par)
    @unpack α, β = par
    return exp(-α-β*N)
end

let 
    Nrange = 0:0.1:10
    RMaturedata = [RMature(N, RickerPar()) for N in Nrange]
    plot(Nrange, RMaturedata)
    xlabel!("N")
    ylabel!("p(Nt)")
end

function RImmature(N, par)
    @unpack α, β = par
    return exp(-α-β*N)
end

let 
    τrange = 0:1:6
    gdata = [calc_g(RickerPar(a=75.0,τ=τval)) for τval in τrange]
    scatter(τrange, gdata, label="")
    xlabel!("τ")
    ylabel!("g(τ)")
end

let 
    τrange = 0:1:6
    pdata = [0.8^(τval+1) for τval in τrange]
    scatter(τrange, pdata, label="")
    xlabel!("τ")
    ylabel!("p(τ)")
end