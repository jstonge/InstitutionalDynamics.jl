using Distributed

@everywhere include("sourcesink2.jl") # import helpers at the same time

# I vs β ------------------------------------------------------------------

@everywhere function get_Iℓ(β, ξ, α, γ, ρ, η, b, c, μ, L, t_max, lvl_1_inf)
  p = [β, ξ, α, γ, ρ, η, b, c, μ]
  sol = run_source_sink2(p, L=L, lvl_1_inf=lvl_1_inf)
  res, res_prop = parse_sol(sol)

  Iℓ = [res[ℓ][end] for ℓ in 1:L]
  Pℓ = [res_prop[ℓ][end] for ℓ in 1:L]
  return Iℓ, Pℓ
end

"""
Run model for many betas here with L=4.
Returns Iℓ (vals) and Pℓ (prop_vals)
"""
function get_base_I_vs_β(βs, L, t_max)
  create_cache()

  cache_I_vs_β = ".cache/res_I&P_vs_β_η_0001_levs.jld"
  if isfile(cache_I_vs_β)
    Iℓ, Pℓ = load(cache_I_vs_β, "I", "P")
  else
    println("P_vs_β_η_0001_levs is not cached yet. We need to run it first.")
    
    @everywhere begin
      L = L
      t_max = t_max
    end

    res_I_vs_β_η_0001_levs = pmap(β -> get_Iℓ(β, 1, 1, 1., 0.05, 0.001, -1., 1., 0.0001, L, t_max, false), βs)  # β, ξ, α, γ, ρ, η, b, c, μ
    Iℓ = [[res_I_vs_β_η_0001_levs[i][1][ℓ] for i in 1:length(βs)] for ℓ in 1:L] # double list comprehension of infected size β x L
    Pℓ = [[res_I_vs_β_η_0001_levs[i][2][ℓ] for i in 1:length(βs)] for ℓ in 1:L]
    save(cache_I_vs_β, "I", Iℓ, "P", Pℓ)
  end
  return Iℓ, Pℓ
end

function run_I₁(βs)
  cache_I1 = ".cache/I1.jld"
  if isfile(cache_I1)
    I₁ = load(".cache/I1.jld", "I1")
  else

    @everywhere begin
      L = 1
      t_max = 10000
    end

    println("I₁ is not cached yet. We need to run it first.")
    I₁ = pmap(β -> get_Iℓ(β, 1, 1, 1., 0.05, 0.005, -1., 1., 0.0001, L, t_max, false), βs)
    I₁ = [x[1][1] for x in I₁] #!TODO: ask Giulio why pmap returns a vector of tuples here?
    save(".cache/I1.jld", "I1", I₁)
  end

  return I₁
end

"""
Make that mapping from ℓ = 1 to higher levels without running things again.
"""
function find_I_ℓ(ℓ::Int, βs::Vector{Float64})
  I₁ = run_I₁(βs) # we always need I₁ for all levels

  if ℓ == 1
    I_ℓ, βsℓ = I₁, βs
  else
    indxs = [findfirst(x -> x >= ℓ*βs[i], βs) for i in 1:length(βs)]
    #     This mapping works ^^^^^^^^  b/c we assume that inst min spread by β_0ℓ^−α
    #     Here we are assuming that α=1, else we would need to multiply ℓ^α
    indxs = indxs[isnothing.(indxs) .== false]
    
    free_zone_ℓ = 1:ℓ:(indxs[1]-1)

    βsℓ = ℓ*βs
    βsℓ = βsℓ[1:findlast(x -> x <= βs[end], βsℓ)]
    βsℓ = [βs[free_zone_ℓ]; βsℓ]

    I_ℓ = zeros(length(indxs))
    I_ℓ = I₁[1:length(indxs)]
    I_ℓ = [zeros(length(free_zone_ℓ)); I_ℓ]
  end

  return βsℓ, I_ℓ
end

#!TODO
# I vs ρ -----------------------------------------------------------------

function get_base_I_vs_ρ(ρs, L, t_max)
  create_cache()

  cache_I_vs_ρ = ".cache/res_I&P_vs_ρ_η_0005_levs.jld"
  if isfile(cache_I_vs_ρ)
    Iℓ, Pℓ = load(cache_I_vs_ρ, "I", "P")
  else
    println("P_vs_ρ_η_0005_levs is not cached yet. We need to run it first.")
    
    @everywhere begin
      L = L
      t_max = t_max
    end

    res_I_vs_ρ_η_0005_levs = pmap(ρ -> get_Iℓ(0.12, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, L, t_max, false), ρs)  # β, ξ, α, γ, ρ, η, b, c, μ

    Iℓ = [[res_I_vs_ρ_η_0005_levs[i][1][ℓ] for i in 1:length(ρs)] for ℓ in 1:L]
    Pℓ = [[res_I_vs_ρ_η_0005_levs[i][2][ℓ] for i in 1:length(ρs)] for ℓ in 1:L]
    save(cache_I_vs_ρ, "I", Iℓ, "P", Pℓ)
  end
  return Iℓ, Pℓ
end

function find_I_ℓ_full(ℓ::Int,  ρs::Vector{Float64})
  # We need to run ρ "fully" bc ρ is independent of institution level, unlike β.
  # There is no way of predicting the effect except by running everything L times. 
  cache_Iℓ_full = ".cache/res_I_vs_ρ_levs$(ℓ).jld"
  if isfile(cache_Iℓ_full)
    I_ℓ = load(cache_Iℓ_full, "I$(ℓ)")
  else
    # ℓ=3
    # t_max = 10_000
    @everywhere begin
      L = 1
      t_max = 10000
    end

    println("res_I_vs_ρ_levs$(ℓ) not cached. Need to run it.")   
    res_I_vs_ρ_lev = pmap(ρ -> get_Iℓ(0.12/ℓ, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, L, t_max, false), ρs)
    I_ℓ = [res_I_vs_ρ_lev[i][1][1] for i in 1:length(ρs)]
    save(cache_Iℓ_full, "I$(ℓ)", I_ℓ)
  end
  return I_ℓ
end

ρs=10 .^ [-3.6:0.02:-0.;]


# scatter(ρs, I_ℓs, xlabel = L"\rho", ylabel = L"\tilde{I}_{\ell}", width = 4., label =:none,
#         palette = palette(:Blues)[[3;5;7;9;]], xscale=:lin)#, xlims = [0.06,0.4])
# savefig("tilde{I}_l_vs_rho.pdf")


# function plot_free_riding_vs_rho_eta_0005(;ρs = 10 .^[-3.6:0.02:-0.;], out="free_riding_vs_rho_eta_0005.pdf")
plot_setup()

L, t_max = 4, 10000 

Iℓs, Pℓ = get_base_I_vs_ρ(ρs, L, t_max)
I_ℓs = [find_I_ℓ_full(ℓ, ρs) for ℓ=1:4]

# scatter(ρs, Iℓs[1], xlabel = L"\rho", ylabel = L"I_\ell", width = 4.,# xlims = [0.06,0.4],
        # labels=:none, legend=:bottomright, palette = palette(:Blues)[[3;5;7;9;]], xscale=:log10);
# plot!(ρs, [sum([Iℓs[ℓ].*Pℓ[ℓ] for ℓ in 1:L])], label=:none, ls=:dash, color=:black, width = 3);
# annotate!(0.015,0.665, text(L"η = 0.005",18))
# savefig("I_l_vs_rho_eta_0005.pdf")

"""
params
======
 * ℓ: chosen ℓ that index global Iℓ
 * I_ℓ: I₁, ... I_L
"""
function get_norm_diff(I_ℓ::Vector{Float64}, Iℓ, ℓ::Int)
  jℓ = findfirst(x -> x > 0.001, I_ℓs) # find index critical point, or where eq prevalence != 0 
  Iℓ_global = Iℓ[1:ℓ:end][(jℓ+1):end]
  tildeI = I_ℓ[(jℓ+1):end]
  repℓ = ℓ == 1 ? 0 : 1
  norm_diff = [repeat([repℓ],jℓ);  (Iℓ_global .- tildeI) ./ (Iℓ_global .+ tildeI)]
  return norm_diff
end

norm_diff_ℓs = [get_norm_diff(I_ℓ, Iℓ, ℓ) for (I_ℓ, Iℓ, ℓ) in zip(I_ℓs, Iℓs, 1:L)]

ρ_cr_l1 = ρs[findfirst(x -> x > 0.001, I_ℓs[1])+1]
scatter(ρs, norm_diff_ℓs, 
        xlabel = L"\rho", palette = palette(:Blues)[[3;5;7;9;]],
        xlims = [ρ_cr_l1, 1.02], ylims = [-1.14,1.14], 
        xscale=:lin, ylabel = L"(I_{\ell} - \tilde{I}_{\ell})/(I_{\ell} + \tilde{I}_{\ell})", 
        legend=:bottomright, width = 4., label =:none);
hline!([0.], ls =:dash, color =:black, width = 2., label=:none);
annotate!(0.5,1.13, text(L"η = 0.005",18));
annotate!(0.5,-0.5, text(L"\textrm{free–riding}",16))
isnothing(out) ? p : savefig(out) 
# end




# prototyping -------------------------------------------------------------------------------

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
# res_I_vs_β_η_0001_levs = pmap(β -> get_Iℓ(β, 1, 1, 1., 0.05, 0.001, -1., 1., 0.0001, L, t_max, false), βs)  # β, ξ, α, γ, ρ, η, b, c, μ
# Iℓ = [[res_I_vs_β_η_0001_levs[i][1][ℓ] for i in 1:length(βs)] for ℓ in 1:L]
# Pℓ = [[res_I_vs_β_η_0001_levs[i][2][ℓ] for i in 1:length(βs)] for ℓ in 1:L]
# save("res_I&P_vs_β_η_0001_levs.jld", "I", Iℓ, "P", Pℓ)
Iℓ, Pℓ = load(".cache/res_I&P_vs_β_η_0001_levs.jld", "I", "P")
scatter(βs, Iℓ, xlabel = L"\beta", ylabel = L"I_\ell", width = 4., xlims = [0.06,0.4],
        labels=:none, legend=:bottomright, palette = palette(:Blues)[[3;5;7;9;]]);
plot!(βs, [sum([Iℓ[ℓ].*Pℓ[ℓ] for ℓ in 1:L])], label=:none, ls=:dot, color=:black, width = 3.);
annotate!(0.23,0.84, text(L"η = 0.05",18))
# savefig("I_l_vs_beta_eta_005.pdf")
# to load back the results from .jld

# This is the script to make that mapping from ℓ = 1 to higher levels without running things again.
# It is messy and not so easy generalizable for arithmetic reasons haha (maybe we can do it later)
# I₁ = pmap(β -> get_Iℓ(β, 1, 1, 1., 0.05, 0.005, -1., 1., 0.0001, 1, t_max, false), βs)
I₁ = load(".cache/res_I_vs_ρ_levs1.jld")
I₁ = I₁["I1"]
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

# save("res_I_vs_β_levs.jld", "I1", I₁, "beta1", βs, "I2", I₂, "beta2", βs2, "I3", I₃, "beta3", βs3, "I4", I₄, "beta4", βs4)
scatter([βs,βs2,βs3,βs4], [I₁,I₂,I₃,I₄], xlabel = L"\beta", ylabel = L"\tilde{I}_{\ell}", width = 4., label =:none,
        palette = palette(:Blues)[[3;5;7;9;]], xlims = [0.06,0.4])
# savefig("tilde{I}_l_vs_beta.pdf")
# to load back the results from .jld
I₁, I₂, I₃, I₄ = load("res_I_vs_ρ_levs.jld", "I1", "I2", "I3", "I4")

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
# savefig("free_riding_vs_beta_eta_005.pdf")







### I vs ρ ###
ρs = 10 .^ [-3.6:0.02:-0.;]
res_I_vs_ρ_η_0005_levs = pmap(ρ -> get_Iℓ(0.12, 1., 1., 1., ρ, 0.005, -1., 1., 0.0001, L, t_max, false), ρs)  # β, ξ, α, γ, ρ, η, b, c, μ
Iℓ = [[res_I_vs_ρ_η_0005_levs[i][1][ℓ] for i in 1:length(ρs)] for ℓ in 1:L]
Pℓ = [[res_I_vs_ρ_η_0005_levs[i][2][ℓ] for i in 1:length(ρs)] for ℓ in 1:L]
# save("res_I&P_vs_ρ_η_0005_levs.jld", "I", Iℓ, "P", Pℓ)
scatter(ρs, Iℓ, xlabel = L"\rho", ylabel = L"I_\ell", width = 4.,# xlims = [0.06,0.4],
        labels=:none, legend=:bottomright, palette = palette(:Blues)[[3;5;7;9;]], xscale=:log10);
plot!(ρs, [sum([Iℓ[ℓ].*Pℓ[ℓ] for ℓ in 1:L])], label=:none, ls=:dash, color=:black, width = 3);
annotate!(0.015,0.665, text(L"η = 0.005",18))
# savefig("I_l_vs_rho_eta_0005.pdf")
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
# savefig("tilde{I}_l_vs_rho.pdf")
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
# savefig("free_riding_vs_rho_eta_0005.pdf")


