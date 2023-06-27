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

"""
params
======
 * ℓ: chosen ℓ that index global Iℓ
 * I_ℓ: I₁, ... I_L
"""
function get_norm_diff(I_ℓ::Vector{Float64}, Iℓ, ℓ::Int)
  jℓ = findfirst(x -> x > 0.001, I_ℓ) # find index critical point, or where eq prevalence != 0 
  Iℓ_global = Iℓ[1:ℓ:end][(jℓ+1):end]
  tildeI = I_ℓ[(jℓ+1):end]
  repℓ = ℓ == 1 ? 0 : 1
  norm_diff = [repeat([repℓ],jℓ);  (Iℓ_global .- tildeI) ./ (Iℓ_global .+ tildeI)]
  return norm_diff
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

# scatter(ρs, Iℓ, xlabel = L"\rho", ylabel = L"I_\ell", width = 4.,# xlims = [0.06,0.4],
#         labels=:none, legend=:bottomright, palette = palette(:Blues)[[3;5;7;9;]], xscale=:log10);
# plot!(ρs, [sum([Iℓ[ℓ].*Pℓ[ℓ] for ℓ in 1:L])], label=:none, ls=:dash, color=:black, width = 3);
# annotate!(0.015,0.665, text(L"η = 0.005",18))
# savefig("I_l_vs_rho_eta_0005.pdf")

# scatter(ρs, I_ℓs, xlabel = L"\rho", ylabel = L"\tilde{I}_{\ell}", width = 4., label =:none,
#         palette = palette(:Blues)[[3;5;7;9;]], xscale=:lin)#, xlims = [0.06,0.4])
# savefig("tilde{I}_l_vs_rho.pdf")


function plot_free_riding_vs_rho_eta_0005(;ρs = 10 .^[-3.6:0.02:-0.;], out="free_riding_vs_rho_eta_0005.pdf")
  # ρs=10 .^ [-3.6:0.02:-0.;]
  plot_setup()

  L, t_max = 4, 10000 
  Iℓs, Pℓ = get_base_I_vs_ρ(ρs, L, t_max)
  I_ℓs = [find_I_ℓ_full(ℓ, ρs) for ℓ=1:4]
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
end