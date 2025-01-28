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
      "-g"
      arg_type = Float64
      default = 1.
      help = "Recovery rate gamma, i.e. rate at which adopters loose behavioral trait"
      "-r"
      arg_type = Float64
      default = 0.1
      help = "Global behavioral diffusion ρ (allows the behaviour to spread between groups)"
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

function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)
  G = zeros(L, n+1)
  for ℓ in 1:L
    for _ in 1:(floor(M/L))
      i = sum(rand(Distributions.Binomial(1, p), n)) # how many total adopters?
      G[ℓ, i+1] += 1                   # everytime combination G[ℓ,i], count +1
    end
  end

  G = G ./ sum(G) # normalized by tot number of groups
  
  # ArrayPartition are nice because we can still access the level such as x[ℓ][i]
  return RecursiveArrayTools.ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
end

function source_sink!(du, u, p, t)
  G, L, n = u, length(u.x), size(u.x[1])[1]
  β, γ, ρ, b, c, μ = p
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
    du.x[ℓ][i] = -γ*n_adopt*G.x[ℓ][i] - (ℓ-1)*β*(n_adopt+R)*(gr_size-n_adopt)*G.x[ℓ][i]

    n_adopt > 0 && ( du.x[ℓ][i] += β*(ℓ-1)*(n_adopt-1+R)*(gr_size-n_adopt+1)*G.x[ℓ][i-1])
    n_adopt < gr_size && ( du.x[ℓ][i] +=  γ*(n_adopt+1)*G.x[ℓ][i+1] )

    # Group selection process
    ℓ > 1 && ( du.x[ℓ][i] += ρ*G.x[ℓ-1][i]*(Z[ℓ] / Z[ℓ-1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ-1] / Z[ℓ]+μ) )
    ℓ < L && ( du.x[ℓ][i] += ρ*G.x[ℓ+1][i]*(Z[ℓ] / Z[ℓ+1] + μ) - ρ*G.x[ℓ][i]*(Z[ℓ+1] / Z[ℓ]+μ) )
  end
end

function run_source_sink(p)
  n, M = 20, 1000
  u₀ = initialize_u0(n=n, L=6, M=M, p=0.05)
  tspan = (1.0, 4000)
  prob = OrdinaryDiffEq.ODEProblem(source_sink!, u₀, tspan, p)
  
  return OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end

function main()
    args = parse_commandline()
    
    if isnothing(args["db"])
      β = args["beta"]
      γ = args["g"]
      ρ = args["r"]
      b = args["b"]
      c = args["c"]
      μ = args["m"]
      
      p = [β, γ, ρ, b, c, μ]
      sol = run_source_sink(p)
      write_sol2txt("$(args["o"])/sourcesink1_$(join(p, "_")).txt", sol) 

    else
      
      db = SQLite.DB(args["db"])
      con = SQLite.DBInterface.execute(db, """SELECT * from sourcesink1 LIMIT $(args["O"]), $(args["L"])""") |> DataFrames.DataFrames.DataFrame
    
      for row in eachrow(con)
        β = row["beta"]
        γ = row["gamma"]
        ρ = row["rho"]
        b = row["b"]
        c = row["cost"]
        μ = row["mu"]
      
        p = [β,  γ, ρ, b, c, μ]
        sol = run_source_sink(p)    
        write_sol2txt("$(args["o"])/sourcesink1_$(join(p, "_")).txt", sol) 
      end
    end  
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end


# prototyping -------------------------------------------------------------------------------

# using Plots

# n, M = 20, 1000
# u₀ = initialize_u0(n=n, L=6, M=M, p=0.001)
# t_max = 1000
# tspan = (0., t_max)

# β, γ, ρ, b, c, μ = 0.22, 0.9, 0.1, 0.18, 1.85, 0.0001
# p  = [β, γ, ρ, b, c, μ]
# prob1 = ODEProblem(source_sink!, u₀, tspan, p)
# sol1 = solve(prob1, DP5(), saveat=1, reltol=1e-8, abstol=1e-8)

# inst_level, inst_level_prop = parse_sol(sol1)  # params: β, γ, ρ, b, c, μ, δ

# function plot_scatter(res::Dict, res_prop::Dict; plot_prop=false)
#   L = length(res)
#   tmax = length(res[L[1]])
#   if plot_prop
#     scatter(1:tmax, [res_prop[i] for i in 1:L], xaxis=:log, legendtitle= "level", 
#           legend=:outertopright, labels = collect(1:L)', palette = palette(:Blues)[3:2:3*(L-1)],
#           markerstrokewidth = 0, markersize = 3.)
#   else 
#     scatter(1:length(res[1]), [res[i] for i in 1:L], xaxis=:log, legendtitle= "level", 
#           legend=:outertopright, labels = collect(1:L)', palette = palette(:Reds)[3:2:3*(L-1)],
#           markerstrokewidth = 0, markersize = 3.)
#     global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:tmax]
#     plot!(1:tmax, global_freq, linestyle=:dash, color=:black, width = 1.5, label = "global") 
#   end
# end

# L = 6
# xs = 1:t_max
# ys = [inst_level_prop[i] for i in 1:L]

# plot_prop
# scatter(xs, ys, xaxis=:log, legendtitle= "level", 
#           legend=:outertopright, labels = collect(1:L)', palette = palette(:Blues),
#           markerstrokewidth = 0, markersize = 3.)

