module InstitutionalDynamics

import StatsBase
import OrdinaryDiffEq
import DataFrames
import SQLite
import Plots

include("sourcesink1.jl")


function sourcesink(;β=0.07, γ=1., ρ=.1, b=.18, c=1.05, μ=.0001)
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


end # module InstitutionalDynamics
