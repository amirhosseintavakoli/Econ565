# Packages

# See the below link for Distribution
# https://juliastats.github.io/Distributions.jl/stable/starting/
using LinearAlgebra
using Random, Distributions
using Optim
using NLsolve
using DataFrames
using DelimitedFiles

# Parameters

J = 2                 # Number of products
N = 100               # Number of individuals
M = 800               # Number of markets
Jchar = 1             # Number of product characteristics
iter = 100            # Number of iterations
β_0 = 5
β_x = 2
# σ_d = 1 # or
σ_d = 3
α = 1
γ_o = 1
γ_x = .5
γ_w = .25
σ_w = .25
σ_c = .25
Y = zeros(2*M, 1)
X1 = zeros(2*M , 3)
β = zeros(3, iter)

## equilibriumPrice_temp ##
# price is an equilibrium variable
# I used FOC for the logit model to compute the price
function equilibriumPrice_temp(price , Variables, ϵ)
    X = Variables[1,:]
    ξ = Variables[2,:]
    w = Variables[3,:]
    ω = Variables[4,:]
    u = [β_0 .+ β_x*X + σ_d*ξ - α*price'; 0] .+ ϵ
    u_max = maximum(u, dims = 1)
    choice_mat = (u_max .<= u)
    s = sum( choice_mat, dims = 2 )/N
    c = exp.(γ_o .+ γ_x*X + σ_c*ξ + γ_w*w + σ_w*ω) # marginal cost (observed)
    equ = price' - c - 1 ./ ( α*(1 .- s[1:2]) )
    return equ
end

# Data Generating Process dgp ##
# This function is used to generate the observed data
# Observed data: p (prices), X(characteristics), c(marginal cost), s(market share), w(input characteristics)
function dgp(i)
    # Variables
    d1 = Normal()
    Random.seed!(200 + i)           # Setting the seed
    Variables = rand(d1, 4, J)
    X = Variables[1,:]
    ξ = Variables[2,:]
    w = Variables[3,:]
    ω = Variables[4,:]
    d2 = Gumbel()
    Random.seed!(300 + i)           # Setting the seed
    ϵ = rand(d2 , J+1 , N)          # error term
    # Find the equilibrium prices
    equilibriumPrice(price) = equilibriumPrice_temp(price, Variables, ϵ ) # Auxiliary function
    initial_p = [.1  .1]
    solution = nlsolve( equilibriumPrice, initial_p)
    p = solution.zero'

    u = [β_0 .+ β_x*X + σ_d*ξ - α*p; 0] .+ ϵ
    u_max = maximum(u, dims = 1)
    choice_mat = (u_max .<= u)
    s = sum( choice_mat, dims = 2)/N
    c = exp.(γ_o .+ γ_x*X + σ_c*ξ + γ_w*w + σ_w*ω)
    return X, w, p, c, s
end

## OLS function ##
function OLS(Y , X)

    #if rank( X'*X ) == 3
        β = ( (X'*X)^(-1) ) * X'*Y
    #end
    #return [0 ; 0 ; 0]
    var_β = (X'*X)^(-1)
    se = σ_d.*sqrt.( diag(var_β) )
    return β , se
end

## MAIN BODY

Y = zeros(2*M*iter, 1)
X1 = zeros(2*M*iter , 5)

for m = 1:M*iter
    try
        X, w, p , c, s = dgp(m)
        y = log.(s[1:2]) .- log( s[3] )
        x1 = [[m;m] X p w c]
        Y[2*m-1:2*m,:] = y
        X1[2*m-1:2*m,:] = x1
        catch err
    end
end

# writedlm( "data_market.csv",  [Y X1], ',') # for σ_d = 1
writedlm( "data_market3.csv",  [Y X1], ',') # for σ_d = 3
