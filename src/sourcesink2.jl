using Distributions, OrdinaryDiffEq, RecursiveArrayTools, ArgParse
include("helpers.jl")

function parse_commandline()
    s = ArgParse.ArgParseSettings()
  
    ArgParse.@add_arg_table! s begin
        "--db"
        help = "Use Database to query parameters"
        "-L"
        arg_type = Int
        default = 5
        help = "LIMIT of rows"
        "-O"
        arg_type = Int
        default = 0
        help = "The OFFSET clause after LIMIT specifies how many rows to skip at the beginning of the result set."
        "--beta"
        arg_type = Float64
        default = 0.07
        help = "Spreading rate from non-adopter to adopter beta"
        "--xi"
        arg_type = Float64
        default = 1.
        help = "Simple-complex contagion parameter"
        "-a"
        arg_type = Float64
        default = 0.5
        help = "Negative benefits alpha"
        "-g"
        arg_type = Float64
        default = 1.
        help = "Recovery rate gamma, i.e. rate at which adopters loose behavioral trait"
        "-r"
        arg_type = Float64
        default = 0.1
        help = "rate between groups to spread the contagion"
        "-e"
        arg_type = Float64
        default = 0.1
        help = "rate between groups to spread the institution level."
        "-b"
        arg_type = Float64
        default = 0.18
        help = "Group benefits b"
        "-c"
        arg_type = Float64
        default = 1.05
        help = "Institutional cost c"
        "-m"
        arg_type = Float64
        default = 1e-4
        help = "Noise u"
        "-o"
        default = "."
        help = "Output file for results"
      end
  
    return ArgParse.parse_args(s)
end

function initial_cond(;n::Int=20, L::Int=4, M::Int=1000000, p::Float64=0.0001, lvl_1::Bool=false)
  G = zeros(L, n+1)

  if lvl_1 # 99% of population at the lowest level
    M /= 10
    for _ in 1:99*M
      i = sum(collect(rand(Binomial(1, p), n))) # how many total infectees?
      G[1, i+1] += 1 # everytime combination [1,i], count +1
    end
    for _ in 1:M
      ℓ = rand(2:L) # pick a level
      i = sum(collect(rand(Binomial(1, p), n))) # how many total infectees?
      G[ℓ, i+1] += 1 # everytime combination [ℓ,i], count +1
    end
    G = G ./ (100*M) # normalized by tot number of groups
    return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L])) 
  else # uniform distribution of population across levels
      for _ in 1:M
        ℓ = rand(1:L) # pick a level
        i = sum(collect(rand(Binomial(1, p), n))) # how many total infectees?
        G[ℓ, i+1] += 1 # everytime combination [ℓ,i], count +1
      end
    G = G ./ M # normalized by tot number of groups
    return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
  end
end

g(x; ξ=1) = x^ξ # function to choose between linear (ξ = 1) and nonlinear contagion (ξ ≠ 1)

@doc raw"""
Key parts of the model:
```
# eq.2 and eq.3
for ℓ in 1:L
  n_infect = collect(0:(n-1))
  R += sum(ρ * n_infect .* G.x[ℓ])
  pop[ℓ] = sum(G.x[ℓ])
  Z[ℓ] = pop[ℓ] > 0 ? sum(exp.(-b*n_infect .- c*(ℓ-1)) .* G.x[ℓ])/pop[ℓ] : 0. 
end

# eq.1 and eq.4
for ℓ = 1:L, i = 1:n
  n_infect, gr_size = i-1, n-1
  # Diffusion
  du.x[ℓ][i] = -γ*n_infect*G.x[ℓ][i] - β*(ℓ^-α)*g(n_infect+R, ξ=ξ)*(gr_size-n_infect)*G.x[ℓ][i]
  n_infect > 0 && ( du.x[ℓ][i] += β*(ℓ^-α)*g(n_infect-1+R, ξ=ξ)*(gr_size-n_infect+1)*G.x[ℓ][i-1])
  n_infect < gr_size && ( du.x[ℓ][i] += γ*(n_infect+1)*G.x[ℓ][i+1] )
  # Selection
  ℓ > 1 && ( du.x[ℓ][i] += η*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - η*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ] + μ) )
  ℓ < L && ( du.x[ℓ][i] += η*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - η*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ] + μ) )
end
```

See the paper for the description.
"""
function dynamics!(du, u, p, t)
  G, L, n = u, size(u.x,1), size(first(u.x),1)
  β, ξ, α, γ, ρ, η, b, c, μ = p
  Z, pop, R = zeros(L), zeros(L), 0.

  # Mean-field coupling and fitness values
  for ℓ in 1:L
    n_infect = collect(0:(n-1))
    R += sum(ρ * n_infect .* G.x[ℓ])
    pop[ℓ] = sum(G.x[ℓ])
    Z[ℓ] = pop[ℓ] > 0 ? sum(exp.(-b*n_infect .- c*(ℓ-1)) .* G.x[ℓ])/pop[ℓ] : 0. 
  end
    
  for ℓ = 1:L, i = 1:n
    n_infect, gr_size = i-1, n-1
    # Diffusion
    du.x[ℓ][i] = -γ*n_infect*G.x[ℓ][i] - β*(ℓ^-α)*g(n_infect+R, ξ=ξ)*(gr_size-n_infect)*G.x[ℓ][i]
    n_infect > 0 && ( du.x[ℓ][i] += β*(ℓ^-α)*g(n_infect-1+R, ξ=ξ)*(gr_size-n_infect+1)*G.x[ℓ][i-1])
    n_infect < gr_size && ( du.x[ℓ][i] += γ*(n_infect+1)*G.x[ℓ][i+1] )
    # Selection
    ℓ > 1 && ( du.x[ℓ][i] += η*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - η*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ] + μ) )
    ℓ < L && ( du.x[ℓ][i] += η*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - η*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ] + μ) )
  end
end

function run_source_sink2(p; L::Int=4, t_max::Int=20000, perc_inf::Float64=p₀, lvl_1::Bool=false)
  n, M = 20, 1000000
  u₀ = initial_cond(n=n, L=L, M=M, p=perc_inf, lvl_1=lvl_1)
  tspan = (1, t_max)
  
  prob = ODEProblem(dynamics!, u₀, tspan, p)
  return solve(prob, Tsit5(), saveat = 1, reltol=1e-8, abstol=1e-8)
end

function main()
  args = parse_commandline()
  
  # Init
  if isnothing(args["db"])
    
    β = args["beta"]
    ξ = args["xi"]
    α = args["a"]
    γ = args["g"]
    ρ = args["r"]
    η = args["e"]
    b = args["b"]
    c = args["c"]
    μ = args["m"]
    
    p = [β, ξ, α, γ, ρ, η, b, c, μ]
    sol = run_source_sink2(p)
    write_sol2txt("$(args["o"])/sourcesink2_$(join(p, "_")).txt", sol) 
  else
    db = SQLite.DB(args["db"])
    con = SQLite.DBInterface.execute(db, """SELECT * from sourcesink2 LIMIT $(args["L"]) OFFSET $(args["O"])""") |> DataFrames.DataFrame

    for row in eachrow(con)
      
      β = row["beta"]
      ξ = row["xi"]
      α = row["alpha"]
      γ = row["gamma"]
      ρ = row["rho"]
      η = row["eta"]
      b = row["b"]
      c = row["cost"]
      μ = row["mu"]
      
      p = [β, ξ, α, γ, ρ, η, b, c, μ]
      sol = run_source_sink2(p)
      write_sol2txt("$(args["o"])/sourcesink2_$(join(p, "_")).txt", sol)
    end
  end  
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
