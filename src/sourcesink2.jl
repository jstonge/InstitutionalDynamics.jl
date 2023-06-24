import Distributions
import RecursiveArrayTools
import OrdinaryDiffEq
using Distributed
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

@everywhere function run_source_sink2g(p; L=L, perc_inf::Float64=0.001, lvl_1_inf::Bool=false)
  n, M = 20, 1000000
  u₀ = initialize_u0(n=n, L=L, M=M, p=perc_inf, lvl_1_inf=lvl_1_inf)

  tspan = (1, t_max)
  
  # Solve problem
  prob = OrdinaryDiffEq.ODEProblem(source_sink2!, u₀, tspan, p)
  return OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(), saveat = 0.1, reltol=1e-8, abstol=1e-8)
end


@everywhere function get_Iℓ(β, ξ, α, γ, ρ, η, b, c, μ, L, t_max, lvl_1_inf)
  p = [β, ξ, α, γ, ρ, η, b, c, μ]
  sol = run_source_sink2g(p, L=L, lvl_1_inf=lvl_1_inf)
  res, res_prop = parse_sol(sol)

  Iℓ = [res[ℓ][end] for ℓ in 1:L]
  Pℓ = [res_prop[ℓ][end] for ℓ in 1:L]
  return Iℓ, Pℓ
end

function get_fitness_evo(sol)
  L = length(sol.u[1].x)
  n = length(sol.u[1].x[1])
  Z = [zeros(length(sol.t)) for _ in 1:L]
  for i in 1:length(sol.t)
    for ℓ in 1:L
      Z[ℓ][i] = sum(exp.(b*[0:(n-1);] .- c*(ℓ-1)) .* sol.u[i].x[ℓ]) / sum(sol.u[i].x[ℓ])
    end
  end
  return Z
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
using Measures
using Distributed

# Original model 

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


# Institutional free riding

default(legendfont = ("Computer modern", 16),
        tickfont = ("Computer modern", 16),
        guidefontsize = 18, markerstrokewidth=0., markersize = 5,
        linewidth=1, framestyle=:axis,
        titlefontsize=12, grid=:none,
        bottom_margin = 0mm, left_margin = 0mm, right_margin = 0mm)
gr(size=(500,400))

@everywhere begin
  L = 4
  t_max = 10000
end

### I vs β ###
βs = [0.06:0.001:0.4;]
res_I_vs_β_η_0001_levs = pmap(β -> get_Iℓ(β, 1, 1, 1., 0.05, 0.001, -1., 1., 0.0001, L, t_max, false), βs)  # β, ξ, α, γ, ρ, η, b, c, μ
Iℓ = [[res_I_vs_β_η_0001_levs[i][1][ℓ] for i in 1:length(βs)] for ℓ in 1:L]
Pℓ = [[res_I_vs_β_η_0001_levs[i][2][ℓ] for i in 1:length(βs)] for ℓ in 1:L]
save("res_I&P_vs_β_η_0001_levs.jld", "I", Iℓ, "P", Pℓ)
scatter(βs, Iℓ, xlabel = L"\beta", ylabel = L"I_\ell", width = 4., xlims = [0.06,0.4],
        labels=:none, legend=:bottomright, palette = palette(:Blues)[[3;5;7;9;]]);
plot!(βs, [sum([Iℓ[ℓ].*Pℓ[ℓ] for ℓ in 1:L])], label=:none, ls=:dot, color=:black, width = 3.);
annotate!(0.23,0.84, text(L"η = 0.05",18))
savefig("I_l_vs_beta_eta_005.pdf")
# to load back the results from .jld
# Iℓ, Pℓ = load("res_I&P_vs_β_η_005_levs.jld", "I", "P")

# This is the script to make that mapping from ℓ = 1 to higher levels without running things again.
# It is messy and not so easy generalizable for arithmetic reasons haha (maybe we can do it later)
I₁ = pmap(β -> get_Iℓ(β, 1, 1, 1., 0.05, 0.005, -1., 1., 0.0001, 1, t_max, false), βs)
#I₂
indxs = [findfirst(x -> x >= 2*βs[i], βs) for i in 1:length(βs)]
indxs = indxs[isnothing.(indxs) .== false]
βs2 = 2*βs
βs2 = βs2[1:findfirst(x -> x == βs[end], βs2)]
βs2 = [βs[1:2:(indxs[1]-1)];βs2]
I₂ = zeros(length(indxs))
I₂ = I₁[1:length(indxs)]
I₂ = [zeros(length(1:2:(indxs[1]-1)));I₂]
#I₃
indxs = [findfirst(x -> x >= 3*βs[i], βs) for i in 1:length(βs)]
indxs = indxs[isnothing.(indxs) .== false]
βs3 = 3*βs
βs3 = βs3[1:findlast(x -> x <= βs[end], βs3)]
βs3 = [βs[1:3:(indxs[1]-1)];βs3]
I₃ = zeros(length(indxs))
I₃ = I₁[1:length(indxs)]
I₃ = [zeros(length(1:3:(indxs[1]-1)));I₃]
#I₄
indxs = [findfirst(x -> x >= 4*βs[i], βs) for i in 1:length(βs)]
indxs = indxs[isnothing.(indxs) .== false]
βs4 = 4*βs
βs4 = βs4[1:findfirst(x -> x >= βs[end], βs4)]
βs4 = [βs[1:4:(indxs[1]-1)];βs4]
I₄ = zeros(length(indxs))
I₄ = I₁[1:length(indxs)]
I₄ = [zeros(length(1:4:(indxs[1]-1)));I₄]

save("res_I_vs_β_levs.jld", "I1", I₁, "beta1", βs, "I2", I₂, "beta2", βs2, "I3", I₃, "beta3", βs3, "I4", I₄, "beta4", βs4)
scatter([βs,βs2,βs3,βs4], [I₁,I₂,I₃,I₄], xlabel = L"\beta", ylabel = L"\tilde{I}_{\ell}", width = 4., label =:none,
        palette = palette(:Blues)[[3;5;7;9;]], xlims = [0.06,0.4])
savefig("tilde{I}_l_vs_beta.pdf")
# to load back the results from .jld
# I₁, I₂, I₃, I₄ = load("res_I_vs_ρ_levs.jld", "I1", "I2", "I3", "I4")

# computing (Iℓ - Ĩℓ)/(Iℓ + Ĩℓ)
j = findfirst(x -> x > 0.001, I₁)
j2 = findfirst(x -> x > 0.001, I₂)
j3 = findfirst(x -> x > 0.001, I₃)
j4 = findfirst(x -> x > 0.001, I₄)
norm_diff_l1 = [repeat([0],j); (Iℓ[1][(j+1):end] .- I₁[(j+1):end])./(Iℓ[1][(j+1):end] .+ I₁[(j+1):end])]
norm_diff_l2 = [repeat([1],j2); (Iℓ[2][1:2:end][(j2+1):end] .- I₂[(j2+1):end])./(Iℓ[2][1:2:end][(j2+1):end] .+ I₂[(j2+1):end])]
norm_diff_l3 = [repeat([1],j3); (Iℓ[3][1:3:end][(j3+1):end] .- I₃[(j3+1):end])./(Iℓ[3][1:3:end][(j3+1):end] .+ I₃[(j3+1):end])]
norm_diff_l4 = [repeat([1],j4); (Iℓ[4][1:4:end][(j4+1):end] .- I₄[(j4+1):end])./(Iℓ[4][1:4:end][(j4+1):end] .+ I₄[(j4+1):end])]
scatter([βs,βs2,βs3,βs4], [norm_diff_l1,norm_diff_l2,norm_diff_l3,norm_diff_l4], xlabel = L"\beta", palette = palette(:Blues)[[3;5;7;9;]],
        xlims = [βs[j+1],0.4], ylims = [-1.14,1.14],
        ylabel = L"(I_{\ell} - \tilde{I}_{\ell})/(I_{\ell} + \tilde{I}_{\ell})", legend=:bottomright, width = 4., label =:none);
hline!([0.], ls =:dash, color =:black, width = 2., label=:none);
annotate!(0.23,1.13, text(L"η = 0.05",18));
annotate!(0.23,-0.5, text(L"\textrm{free–riding}",16))
savefig("free_riding_vs_beta_eta_005.pdf")


### I vs ρ ###
ρs = 10 .^ [-3.6:0.02:-0.;]
res_I_vs_ρ_η_0005_levs = pmap(ρ -> get_Iℓ(0.12, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, L, t_max, false), ρs)  # β, ξ, α, γ, ρ, η, b, c, μ
Iℓ = [[res_I_vs_ρ_η_0005_levs[i][1][ℓ] for i in 1:length(ρs)] for ℓ in 1:L]
Pℓ = [[res_I_vs_ρ_η_0005_levs[i][2][ℓ] for i in 1:length(ρs)] for ℓ in 1:L]
save("res_I&P_vs_ρ_η_0005_levs.jld", "I", Iℓ, "P", Pℓ)
scatter(ρs, Iℓ, xlabel = L"\rho", ylabel = L"I_\ell", width = 4.,# xlims = [0.06,0.4],
        labels=:none, legend=:bottomright, palette = palette(:Blues)[[3;5;7;9;]], xscale=:log10);
plot!(ρs, [sum([Iℓ[ℓ].*Pℓ[ℓ] for ℓ in 1:L])], label=:none, ls=:dash, color=:black, width = 3);
annotate!(0.015,0.665, text(L"η = 0.005",18))
savefig("I_l_vs_rho_eta_0005.pdf")
# to load back the results from .jld
# Iℓ, Pℓ = load("res_I&P_vs_ρ_η_0005_levs.jld", "I", "P")

res_I_vs_ρ_lev1 = pmap(ρ -> get_Iℓ(0.12, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, 1, t_max, false), ρs)
res_I_vs_ρ_lev2 = pmap(ρ -> get_Iℓ(0.12/2, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, 1, t_max, false), ρs)
res_I_vs_ρ_lev3 = pmap(ρ -> get_Iℓ(0.12/3, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, 1, t_max, false), ρs)
res_I_vs_ρ_lev4 = pmap(ρ -> get_Iℓ(0.12/4, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, 1, t_max, false), ρs)
save("res_I_vs_ρ_levs.jld", "I1", I₁, "I2", I₂, "I3", I₃, "I4", I₄)
I₁ = [res_I_vs_ρ_lev1[i][1][1] for i in 1:length(ρs)]
I₂ = [res_I_vs_ρ_lev2[i][1][1] for i in 1:length(ρs)]
I₃ = [res_I_vs_ρ_lev3[i][1][1] for i in 1:length(ρs)]
I₄ = [res_I_vs_ρ_lev4[i][1][1] for i in 1:length(ρs)]
scatter(ρs, [I₁,I₂,I₃,I₄], xlabel = L"\rho", ylabel = L"\tilde{I}_{\ell}", width = 4., label =:none,
        palette = palette(:Blues)[[3;5;7;9;]], xscale=:lin)#, xlims = [0.06,0.4])
savefig("tilde{I}_l_vs_rho.pdf")
# to load back the results from .jld
# I₁, I₂, I₃, I₄ = load("res_I_vs_ρ_levs.jld", "I1", "I2", "I3", "I4")

# computing (Iℓ - Ĩℓ)/(Iℓ + Ĩℓ)
j = findfirst(x -> x > 0.001, I₁)
j2 = findfirst(x -> x > 0.001, I₂)
j3 = findfirst(x -> x > 0.001, I₃)
j4 = findfirst(x -> x > 0.001, I₄)
ρ_cr_l1 = ρs[j+1]
norm_diff_l1 = [repeat([0],j); (Iℓ[1][(j+1):end] .- I₁[(j+1):end])./(Iℓ[1][(j+1):end] .+ I₁[(j+1):end])]
norm_diff_l2 = [repeat([1],j2); (Iℓ[2][(j2+1):end] .- I₂[(j2+1):end])./(Iℓ[2][(j2+1):end] .+ I₂[(j2+1):end])]
norm_diff_l3 = [repeat([1],j3); (Iℓ[3][(j3+1):end] .- I₃[(j3+1):end])./(Iℓ[3][(j3+1):end] .+ I₃[(j3+1):end])]
norm_diff_l4 = [repeat([1],j4); (Iℓ[4][(j4+1):end] .- I₄[(j4+1):end])./(Iℓ[4][(j4+1):end] .+ I₄[(j4+1):end])]
scatter(ρs, [norm_diff_l1,norm_diff_l2,norm_diff_l3,norm_diff_l4], xlabel = L"\rho", palette = palette(:Blues)[[3;5;7;9;]],
        xlims = [ρ_cr_l1,1.02], ylims = [-1.14,1.14], 
        xscale=:lin, ylabel = L"(I_{\ell} - \tilde{I}_{\ell})/(I_{\ell} + \tilde{I}_{\ell})", legend=:bottomright, width = 4., label =:none);
hline!([0.], ls =:dash, color =:black, width = 2., label=:none);
annotate!(0.5,1.13, text(L"η = 0.005",18));
annotate!(0.5,-0.5, text(L"\textrm{free–riding}",16))
savefig("free_riding_vs_rho_eta_0005.pdf")