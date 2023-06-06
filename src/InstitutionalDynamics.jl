module InstitutionalDynamics

import StatsBase
import OrdinaryDiffEq
import DataFrames
import SQLite
import Plots

include("sourcesink1.jl")
include("sourcesink2.jl")
include("sourcesink3.jl")


function sourcesink_m(;β=0.07, γ=1., ρ=.1, b=.18, c=1.05, μ=.0001)
  p = [β,  γ, ρ, b, c, μ]
  sol = run_source_sink(p)
  inst_level, inst_level_prop = parse_sol(sol)
  
  L = 6
  t_max = 1000
  xs = 1:t_max
  ys = [inst_level_prop[i] for i in 1:L]  
  
  return Plots.scatter(xs, ys, xaxis=:log, legendtitle= "level", 
                       legend=:outertopright, labels = collect(1:L)', palette = Plots.palette(:Blues),
                       markerstrokewidth = 0, markersize = 3.)
end

function contagion_m(;β=0.1, ξ=1., α=1., γ=1., ρ=0.1, η=0.05, b=-1., c=1., μ=.0001)
  p = [β, ξ, α, γ, ρ, η, b, c, μ]
  sol = run_source_sink2(p)
  inst_level, inst_level_prop = parse_sol(sol)
  
  L = 4
  t_max = 1000
  xs = 1:t_max
  ys = [inst_level_prop[i] for i in 1:L]  
  
  return Plots.scatter(xs, ys, xaxis=:log, legendtitle= "level", 
                       legend=:outertopright, labels = collect(1:L)', 
                       palette = Plots.palette(:Blues),
                       markerstrokewidth = 0, markersize = 3.)
end

function ressource_m(;β=0.07, γ=.1, ρ=.5, b=.18, c=1.05, μ=.2,  δ=1.0, α=0.2)
  p = [β, γ, ρ, b, c, μ, δ, α]
  sol = run_source_sink3(p)
  inst_level, inst_level_prop = parse_sol(sol)
  
  L = 4
  t_max = 1000
  xs = 1:t_max
  ys = [inst_level_prop[i] for i in 1:L]  
  
  return Plots.scatter(xs, ys, xaxis=:log, legendtitle= "level", 
                       legend=:outertopright, labels = collect(1:L)', palette = Plots.palette(:Blues),
                       markerstrokewidth = 0, markersize = 3.)
end


end # module InstitutionalDynamics
