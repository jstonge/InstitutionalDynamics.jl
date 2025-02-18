using Distributions, RecursiveArrayTools, OrdinaryDiffEq, ArgParse, SQLite, DataFrames, Plots, ColorSchemes

include("helpers.jl")


"""
    parse_commandline()

Function for the commandline interface.
"""
function parse_commandline()
  s = ArgParseSettings()

  #!TODO: Update the params and defaults as needed.
  @add_arg_table! s begin
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
      help = "Imitation rate from non-adopter to adopter β"
      "-g"
      arg_type = Float64
      default = 1.
      help = "Imitation rate from adopter to non-adopter γ"
      "-r"
      arg_type = Float64
      default = 0.1
      help = "Global behavioral imitation ρ (allows a behaviour to spread between groups)"
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
      help = "Endogenous rate of institutional change μ"
      "-a"
      arg_type = Float64
      default = 0.2
      help = "Endogenous rate of individual change α"
      "-o"
      default = "."
      help = "Output file for results"
    end

   return parse_args(s)
end

"""
  initialize_u0(;n, L, M, p)

Function to initialize the model.
"""
function initialize_u0(;n::Int=20, L::Int=6, M::Int=20, p::Float64=0.01)
  G = zeros(L, n+1)
  for _ in 1:M
    ℓ = rand(1:L)
    i = sum(rand(Binomial(1, p), n))
    G[ℓ, i+1] += 1
  end

  G = G ./ M

  return ArrayPartition(Tuple([G[ℓ,:] for ℓ=1:L]))
end

f(x; a1=3) = 1 / (1 + exp(-a1*x)) # individual incentive to cooperate given the (perceived) quality of the public good (∝ ℓ)
s(x; a2=0.3) = a2 > 0 ? (1 - exp(-a2*x)) / (1 - exp(-a2)) : x # dependency of the (perceived) quality of the public good on institutional level
h(x; a3=0.0) = exp(-a3*x) # rescaling to have a (decreasing) level-dependent benefit

# Note 1: Since the difference between the payoff of a defector (D) and of a cooperator (C) 
#         is the same –one unit of "capital"– whatever the group size (n), its composition (i) and the quality of the public good (∝ ℓ);
#         we conflate the imitation due to strategic learning (generating from individual return maximization) 
#         and the imitation due to pure peer-pressure modelling them via a unique constant rate, 
#         one for each strategic change (β: D –> C; γ: C –> D).
# Note 2: In the following, we rescaled time by the rate associated to the institution-dependent individual incentive to cooperate f(...).
# Note 3: The functions g (cost-benefits for groups) and g̃ (fitness function) are taken equal to function f. 
#         The three have similar properties.

function source_sink3!(du, u, p, t)
  G, L, n = u, length(u.x), length(u.x[1])
  β, γ, ρ, b, c, μ, δ, α = p # δ = 1 (δ = 0): (no) resource requirement to upgrade institution
  Z, pop, R = zeros(L), zeros(L), 0.

  # Calculate mean-field coupling and observed fitness landscape
    for ℓ in 1:L
      n_adopt = collect(0:(n-1))
      Z[ℓ]    = sum(f.(b*h(ℓ-1)*n_adopt .- c*(ℓ-1)) .* G.x[ℓ])
      pop[ℓ]  = sum(G.x[ℓ])
      R      += sum(n_adopt .* G.x[ℓ]) # Global diffusion
      pop[ℓ] > 0.0 && ( Z[ℓ] /= pop[ℓ] )
    end

    for ℓ = 1:L, i = 1:n
      # we separate indices from the meaning of the variables, 
      # aka current_level (0,1,2,3) or #adopters vs. accessing Z[ℓ-1] or G.x[ℓ][i-1]
      current_level, n_adopt, gsize = ℓ-1, i-1, n-1
      # Individual selection process
      du.x[ℓ][i] = - α*n_adopt*f(1-s(current_level))*G.x[ℓ][i] - α*(gsize-n_adopt)*f(s(current_level)-1)*G.x[ℓ][i]
      du.x[ℓ][i] += - n_adopt*(gsize-n_adopt)*(β+γ)*G.x[ℓ][i] - ρ*(gsize-n_adopt)*β*R*G.x[ℓ][i] - ρ*n_adopt*γ*(gsize-R)*G.x[ℓ][i]
      n_adopt > 0 && ( du.x[ℓ][i] += α*(gsize-n_adopt+1)*f(s(current_level)-1)*G.x[ℓ][i-1] + β*(n_adopt-1+ρ*R)*(gsize-n_adopt+1)*G.x[ℓ][i-1] )
      n_adopt < gsize && ( du.x[ℓ][i] += α*(n_adopt+1)*f(1-s(current_level))*G.x[ℓ][i+1] + γ*(gsize-n_adopt-1+ρ*(gsize-R))*(n_adopt+1)*G.x[ℓ][i+1] )
      # Group selection process
      ℓ > 1 && ( 
        du.x[ℓ][i] += (f(b*h(current_level)*n_adopt - current_level*c)^δ)*(μ + ρ*Z[ℓ]/Z[ℓ-1])*G.x[ℓ-1][i]
                          - (μ*(f(current_level*c - b*h(current_level)*n_adopt)^δ) + ρ*Z[ℓ-1]/Z[ℓ])*G.x[ℓ][i] 
            )
      ℓ < L && ( 
        du.x[ℓ][i] += (μ*(f(c*ℓ-b*h(ℓ)*n_adopt)^δ) + ρ*Z[ℓ]/Z[ℓ+1])*G.x[ℓ+1][i]
                            - (f(b*h(ℓ)*n_adopt-c*ℓ)^δ)*(μ+ρ*Z[ℓ+1]/Z[ℓ])*G.x[ℓ][i] 
            )
    end
end

function run_source_sink3(p)
  n, M = 20, 1000
  u₀ = initialize_u0(n=n, L=5, M=M, p=0.01)
  t_max = 2000
  tspan = (1.0, t_max)

  # Solve problem
  prob = ODEProblem(source_sink3!, u₀, tspan, p)
  return solve(prob, DP5(), saveat = 1., reltol=1e-8, abstol=1e-8)
end

function main()
  # β, γ, ρ, b, c, μ, δ, α  = 0.07, 0.5, 1, 0.1, 0.18, 1.05, 0.2, 1, 0.2
  args = parse_commandline()

  modelname = "sourcesink3"

  # Init
  if isnothing(args["db"])
    β = args["beta"]
    γ = args["g"]
    ρ = args["r"]
    b = args["b"]
    c = args["c"]
    μ = args["m"]
    δ = args["d"]
    α = args["a"]

    p = [β, γ, ρ, b, c, μ, δ, α]
    println(p)
    sol = run_source_sink3(p)
    write_sol2txt("$(args["o"])$(modelname)_$(join(p, "_")).txt", sol)
  else
    db = SQLite.DB(args["db"])
    con = SQLite.DBInterface.execute(db, """SELECT * from $(modelname) LIMIT $(args["L"]) OFFSET $(args["O"])""") |> DataFrames.DataFrame

    for row in eachrow(con)

      β = row["beta"]
      γ = row["gamma"]
      ρ = row["rho"]
      b = row["b"]
      c = row["cost"]
      μ = row["mu"]
      δ = row["delta"]
      α = row["alpha"]

      p = [β, γ, ρ, b, c, μ, δ, α]
      sol = run_source_sink3(p)
      write_sol2txt("$(args["o"])/$(modelname)_$(join(p, "_")).txt", sol)
    end
  end
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end



# prototyping -------------------------------------------------------------------------------

# # note 1: b/c should be chosen to allow the highest level (hence also all the others) to be worth it at high enough adoption
# #      --> lower bound for b/c (given n = 20 and max{level} = 5): 20*b - 5*c > 0 ⇒ b/c > 0.25
# # note 2: β/γ and ρ/μ decisive for dominance relations between levels
# #      --> heatmap β/γ VS ρ/μ for fixed b/c?

# n, M = 20, 1000
# u₀ = initialize_u0(n=n, L=4, M=M, p=0.01)
# t_max = 1000
# tspan = (0., t_max)

# h(x; a3=0.) = exp(-a3*x)
# # β, γ, ρ, b, c, μ, α = 0.02, 0.02, 0.02, 0.16, 1., 0.1, 0.2
# β, γ, ρ, b, c, μ, α, δ = 0.2, 0.05, 0., 0.2, 1., 0.06, 0.2, 1
# p  = [β, γ, ρ, b, c, μ, δ, α]
# prob = ODEProblem(source_sink3!, u₀, tspan, p)
# sol = solve(prob, Tsit5(), saveat=1, reltol=1e-12, abstol=1e-12)

# δ = 0
# p  = [β, γ, ρ, b, c, μ, δ, α]
# prob1 = ODEProblem(source_sink3!, u₀, tspan, p)
# sol1 = solve(prob1, DP5(), saveat=1, reltol=1e-8, abstol=1e-8)

# inst_level, inst_level_prop = parse_sol(sol)  # params: β, γ, ρ, b, c, μ, δ

# # temporal evolution
# function plot_scatter(res::Dict, res_prop::Dict; plot_prop=false)
#   L = length(res)
#   tmax = length(res[L[1]])
#   if plot_prop
#     scatter(1:tmax, [res_prop[i] for i in 1:L], xaxis=:log, legendtitle= "level", 
#           legend=:outertopright, labels = collect(1:L)', palette = colorschemes[:Blues][3:9],
#           markerstrokewidth = 0, markersize = 3.)
#   else 
#     scatter(1:length(res[1]), [res[i] for i in 1:L], xaxis=:log, legendtitle= "level", 
#           legend=:outertopright, labels = collect(1:L)', palette = colorschemes[:Reds][3:9],
#           markerstrokewidth = 0, markersize = 3.)
#     global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:tmax]
#     plot!(1:tmax, global_freq, linestyle=:dash, color=:black, width = 1.5, label = "global") 
#   end
# end

# # plot results
# plot(
#   plot_scatter(inst_level, inst_level_prop),
#   plot_scatter(inst_level, inst_level_prop, plot_prop = true), 
#   layout=(2,1), size=(800,800))

# # plotting s()
# plot(ℓ -> s(ℓ), 0, 3,  label="a2=0.3")
# plot(ℓ -> s(ℓ, a2=0.1), 0, 3, label="a2=0.1")
# plot!(ℓ -> s(ℓ, a2=0.5), 0, 3, label="a2=0.5")

# # plotting f()
# plot(ℓ -> 1 / (1 + exp(-3*(1-s(ℓ, a2=0.1)))), 0, 3, label="a1=3")
# plot!(ℓ -> 1 / (1 + exp(-3*(1-s(ℓ, a2=0.5)))), 0, 3, label="a1=3")
# ylabel!("cost")

# # plotting f()
# plot(ℓ -> f(1 - s(ℓ)), 0, 3, label="a1=3")
# plot!(ℓ -> f(1 - s(ℓ), a1=0), 0, 3, label="a1=0", ls=:dash)
# plot!(ℓ -> f(1 - s(ℓ), a1=1), 0, 3, label="a1=1")
# plot!(ℓ -> f(1 - s(ℓ), a1=5), 0, 3, label="a1=5")
# title!("Coop given (perceived) quality of the\nPG (a2 fixed at 0.3)")
# ylabel!("low   ⇐    cost    ⇒    high")
# xlabel!("Institutional strength")
# ylims!(0,1)
# annotate!(2.5, 0.55, "line of indifference", color = :blue)


# plot(ℓ -> f(1 - s(ℓ, a2=0.1)), 0, 3, label="a1=3")
# plot!(ℓ -> f(1 - s(ℓ, a2=0.1), a1=0), 0, 3, label="a1=0", ls=:dash)
# plot!(ℓ -> f(1 - s(ℓ, a2=0.1), a1=1), 0, 3, label="a1=1")
# plot!(ℓ -> f(1 - s(ℓ, a2=0.1), a1=5), 0, 3, label="a1=5")
# title!("Coop given (perceived) quality of the\nPG (a2 fixed at 0.1)")
# ylabel!("low   ⇐    cost    ⇒    high")
# xlabel!("Institutional strength")
# ylims!(0,1)
# annotate!(2.5, 0.55, "line of indifference", color = :blue)


###########################
# test reading from sweep #
###########################
#
# fsol = CSV.read("sourcesink3_0.04_0.04_0.09_0.26_1.0_0.1_1.0_0.15.txt", DataFrame; header=["timestep", "L", "value"], 
#                    types=Dict(:timestep => Int, :L => Int, :value => Float32))
# gd = groupby(fsol, [:timestep, :L])
# n = nrow(gd[1])
    
# process functions
# processing1(x, n) = round(sum((collect(0:(n-1)) / n) .* x) / sum(x), digits=7)

# df_agg = combine(gd, :value => (x -> round(sum(x), digits=7)) => :value_prop, 
#                      :value => (x -> iszero(sum(x)) ? 0.0 : processing1(x,n)) => :value)

# scatter(1:maximum(df_agg.timestep), [filter(:L => x -> x ==  i, df_agg).value_prop for i in 1:L], xaxis=:log, legendtitle= "level", 
#           legend=:outertopright, labels = collect(1:L)', palette = colorschemes[:Blues][3:9],
#           markerstrokewidth = 0, markersize = 3.)
#
###########################
