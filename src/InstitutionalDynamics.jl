module InstitutionalDynamics

import StatsBase
import OrdinaryDiffEq
import DataFrames
import SQLite
import Plots

include("sourcesink1.jl")
include("sourcesink2.jl")
include("sourcesink3.jl")
include("p_free_riding.jl")


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


function plot_I_l_vs_beta_eta_005(;βs=[0.06:0.001:0.4;], out="I_l_vs_beta_eta_005.pdf")
  L, t_max = 4, 10_000
  Iℓ, Pℓ = get_base_I_vs_β(βs, L, t_max)
  
  p = scatter(βs, Iℓ, xlabel = L"\beta", 
          ylabel = L"I_\ell", 
          width = 4., xlims = [0.06,0.4],
          labels = :none, 
          legend = :bottomright, 
          palette = palette(:Blues)[[3;5;7;9;]]);
  ys =  [sum([Iℓ[ℓ].*Pℓ[ℓ] for ℓ in 1:L])]
  plot!(βs, ys, label = :none, ls = :dot, color = :black, width = 3.)
  annotate!(0.23,0.84, text(L"η = 0.05",18))
  isnothing(out) ? p : savefig(out)
end


function plot_tildeI_l_vs_beta(;βs=[0.06:0.001:0.4;], out="tildeI_l_vs_beta.pdf")
  L = 4
  βsℓs, I_ℓs = unzip([find_I_ℓ(ℓ, βs) for ℓ=1:L])
  p = scatter(βsℓs, I_ℓs, 
          xlabel = L"\beta", 
          ylabel = L"\tilde{I}_{\ell}", width = 4., label =:none,
          palette = palette(:Blues)[[3;5;7;9;]], xlims = [0.06,0.4])
  
  isnothing(out) ? p : savefig(out)
end


function plot_free_riding(;βs=[0.06:0.001:0.4;], out="free_riding_vs_beta_eta_005.pdf")

  plot_setup()

  L, t_max = 4, 10000 
  Iℓs, Pℓ = get_base_I_vs_β(βs, L, t_max)
  βsℓs, I_ℓs = unzip([find_I_ℓ(ℓ, βs) for ℓ=1:L])
  norm_diff_ℓs = [get_norm_diff(I_ℓ, Iℓ, ℓ) for (I_ℓ, Iℓ, ℓ) in zip(I_ℓs, Iℓs, 1:L)]
  j = findfirst(x -> x > 0.001,  run_I₁(βs))

  p = scatter(βsℓs, norm_diff_ℓs, 
          xlabel = L"\beta", palette = palette(:Blues)[[3;5;7;9;]],
          xlims = [βs[j+1],0.4], ylims = [-1.14,1.14],
          ylabel = L"(I_{\ell} - \tilde{I}_{\ell})/(I_{\ell} + \tilde{I}_{\ell})", 
          legend=:bottomright, width = 4., label =:none);
  
  hline!([0.], ls = :dash, color =:black, width = 2., label=:none);
  annotate!(0.23,1.13, text(L"η = 0.05",18));
  annotate!(0.23,-0.5, text(L"\textrm{free–riding}",16))
  isnothing(out) ? p : savefig(out)  
end


end # module InstitutionalDynamics
