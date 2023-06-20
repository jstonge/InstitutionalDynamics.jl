import Distributions
import RecursiveArrayTools
import OrdinaryDiffEq
import ArgParse
import SQLite
import DataFrames

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

function initialize_u0(;n::Int=20, L::Int=6, M::Int=100, p::Float64=0.001,
  lvl_1_inf::Bool=false)
  G = zeros(L, n+1)

  if lvl_1_inf # 99% of population at the lowest level
    M /= 10
    for _ in 1:99*M
      i = sum(collect(rand(Distributions.Binomial(1, p), n))) # how many total adopters?
      G[1, i+1] += 1 # everytime combination [1,i], count +1
    end
    for _ in 1:M
      ℓ = rand(2:L) # pick a level
      i = sum(collect(rand(Distributions.Binomial(1, 0.01*p), n))) # how many total adopters?
      G[ℓ, i+1] += 1 # everytime combination [ℓ,i], count +1
    end
    G = G ./ (100*M) # normalized by tot number of groups
    return RecursiveArrayTools.ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
  else 
      for _ in 1:M
        ℓ = rand(1:L) # pick a level
        i = sum(collect(rand(Distributions.Binomial(1, p), n))) # how many total adopters?
        G[ℓ, i+1] += 1 # everytime combination [ℓ,i], count +1
      end
    G = G ./ M # normalized by tot number of groups
    return RecursiveArrayTools.ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
  end
end

# function to switch between simple (ξ = 1) and complex contagion (ξ ≠ 1)
g(x; ξ=1.) = x^ξ # for x ∈ ℕ ⇒ ξ = 1: linear growth; 0 < ξ < 1: sublinear growth; ξ > 1: superlinear growth

function source_sink2!(du, u, p, t)
    G, L, n = u, length(u.x), length(first(u.x))
    β, ξ, α, γ, ρ, η, b, c, μ = p
    Z, pop, R = zeros(L), zeros(L), 0.

    # Calculate mean-field coupling and observed fitness landscape
    for ℓ in 1:L
        n_adopt = collect(0:(n-1))
        Z[ℓ]    = sum(exp.(b*n_adopt .- c*(ℓ-1)) .* G.x[ℓ]) 
        pop[ℓ]  = sum(G.x[ℓ])
        R      += sum(ρ * n_adopt .* G.x[ℓ]) 
        pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] ) 
      end
      
      for ℓ = 1:L, i = 1:n
        n_adopt, gr_size = i-1, n-1
        # Diffusion events
        du.x[ℓ][i] = -γ*n_adopt*G.x[ℓ][i] - β*(ℓ^-α)*g(n_adopt+R, ξ=ξ)*(gr_size-n_adopt)*G.x[ℓ][i]
        n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ^-α)*g(n_adopt-1+R, ξ=ξ)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
        n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )
        # Group selection process
        ℓ > 1 && ( du.x[ℓ][i] += η*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - η*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ] + μ) )
        ℓ < L && ( du.x[ℓ][i] += η*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - η*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ] + μ) )
      end
end

function run_source_sink2(p; perc_inf::Float64=0.001, lvl_1_inf::Bool=false)
  n, M = 20, 10000
  L = 4
  u₀ = initialize_u0(n=n, L=L, M=M, p=perc_inf, lvl_1_inf=lvl_1_inf)

  tspan = (1.0, 10000)
  
  # Solve problem
  prob = OrdinaryDiffEq.ODEProblem(source_sink2!, u₀, tspan, p)
  return OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
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

# prototyping -------------------------------------------------------------------------------

using LaTeXStrings
using Plots

params_name = "β", "ξ", "α", "γ", "ρ", "η", "b", "c", "μ"

lvl_1_inf = false
perc_inf = 0.001
p = [0.1, 1., 1., 1., 0.1, 0.05, -1, 1., 0.0001]  # β, ξ, α, γ, ρ, η, b, c, μ
sol = run_source_sink2(p, perc_inf = perc_inf, lvl_1_inf=lvl_1_inf)
res, res_prop = parse_sol(sol)
t_max = 9999

function plot_value(res)
  L = length(res)
  xval = [res[l][1:t_max] for l=1:L]
  plot(xval, 
    xscale=:log, 
    # ylabel = L"\textrm{prevalence}", 
    labels = " " .* string.([1:L;]'),
    width = 3., 
    legendtitle = L"\textrm{level}",
    palette = palette(:Reds)[3:9], 
    legend=:left,
    xticks = 10 .^ [0,1,2,3,4]
  );
  
  global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:t_max]
  
  plot!(1:t_max, global_freq[1:t_max], width = 3,
        color =:black, ls =:dash, label = L"\textrm{global}",
        title = join([params_name[i] * "=" * string.(p)[i] for i in 1:length(params_name)], "  "))
  
end

plot_value(res)

function plot_value_prop(res_prop)
  L = length(res)
  xval = [res_prop[l][1:t_max] for l=1:L]
  plot(xval, xscale=:log, ylabel = L"\textrm{level\ proportion}",
       labels = " " .* string.([1:L;]'),
       width = 3., 
       legendtitle = L"\textrm{level}", 
       palette = palette(:Blues)[3:9], 
       legend=:outerright,
       title = join([params_name[i] * "=" * string.(p)[i] for i in 1:length(params_name)], "  "),
       xticks = 10 .^ [0,1,2,3,4])
end

plot_value_prop(res_prop)