using SymPy

# Define symbolic variables
@syms b D K C α β a b 
@syms τ


g(τ)=a-b*exp(-K*(τ+1))
a0(τ)=D*(α*((1+D)^(τ+1)- g(τ)) - g(τ))
a1(τ)=-((α*C + β*D)*g(τ)) + (1+D)^(τ+1)*(β*D + α*C*g(τ))
a2(τ)=((1+D)^(τ+1)-1)*β*C*g(τ)
x(τ)=(-a1(τ) + sqrt(a1(τ)^2 - 4*a0(τ)*a2(τ))) / (2*a2(τ))
q(τ)=b*D*K - g(τ)*(D + C*g(τ)*x(τ)*exp(-K*(1 + τ)))*log(1+D)

solve(q(τ), τ)

using Plots

function ricker(x, D, C)
    return x * exp(-D-C*x)    
end

function adult(x, α, β)
    return 1-exp(-α-β*x)
end

r1=ricker(1/0.3, 0.1, 0.3)

r2=ricker(r1,0.1, 0.3)

function partial(x, D, C)
    return exp(-D-C*x) * (1 - C*x)
end

function gptilde(x, D, C)
    return x * exp(-D-C*x) + D * (1 - exp(-C*x))
end

let 
    C=0.3
    D=0.1
    α=0.1
    β=2.0
    xrange = 0:0.1:10
    data=[ricker(0.5*x, D, C)/x for x in xrange]
    data2=[adult(x, α, β) for x in xrange]
    data3=[ricker(ricker(0.5*x, D, C), D, C)/x for x in xrange]
    data4=[ricker(0.5*x, D, C) for x in xrange]

    plot(xrange, data, label="Ricker", xlabel="x", ylabel="f(x)")
    plot!(xrange, data2, label="Adult", xlabel="x", ylabel="f(x)")
    plot!(xrange, data3, label="Ricker (2nd)", xlabel="x", ylabel="f(x)")
    plot!(xrange, data4, label="Ricker (0.5)", xlabel="x", ylabel="f(x)")
end

let 
    C=0.3
    D=0.1
    xrange = 0:0.1:10
    data=[ricker(ricker(x, D, C), D, C) for x in xrange]
    plot(xrange, data, xlabel="x", ylabel="f(x)", title="Ricker Function")
end


let 
    C=0.3
    D=0.1
    xrange = 0:0.1:100
    data=[partial(x, D, C) for x in xrange]
    plot(xrange, data, xlabel="x", ylabel="f'(x)")
end

let
    C=0.3
    D=0.1
    xrange = 0:0.1:100
    data=[gptilde(x, D, C) for x in xrange]
    plot(xrange, data, xlabel="x", ylabel="g~(x)")
end