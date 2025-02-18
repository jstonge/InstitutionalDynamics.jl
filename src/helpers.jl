using Printf

function parse_sol(s)
  L = length(s.u[1].x)
  tmax = length(s)-1
  inst_level = Dict()
  inst_level_prop = Dict()
  for ℓ=1:L
    values = []
    values_prop = []
    for t=1:tmax
      n = length(s.u[t].x[ℓ])
      x = s.u[t].x[ℓ]
      out = sum((collect(0:(n-1)) / n) .* x) / sum(x) 
      push!(values, out)
      out = sum(x)
      push!(values_prop, out)
    end
    inst_level[ℓ] = values
    inst_level_prop[ℓ] = values_prop
  end
  return inst_level, inst_level_prop
end

"""
write_sol2txt(path, sol)

Function to write solution to textfile. Nothing to do here.
"""
function write_sol2txt(path, sol)
  L = length(sol.u[1].x)
  open(path, "w") do io  # Use "w" to overwrite instead of appending
      for t in 1:length(sol.u)
          for ℓ in 1:L
              for val in sol.u[t].x[ℓ]
                if val < 0 || val > 1
                  println("$(path), $(val)")
                end
                @printf(io, "%d %d %.14f\n", t, ℓ, val)
              end
          end
      end
  end
end

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))
create_cache() = isdir(".cache") ? nothing : mkdir(".cache")

function plot_setup()
  default(legendfont = ("Computer modern", 16),
          tickfont = ("Computer modern", 16),
          guidefontsize = 18, markerstrokewidth=0., markersize = 5,
          linewidth=1, framestyle=:axis,
          titlefontsize=12, grid=:none,
          bottom_margin = 0mm, left_margin = 0mm, right_margin = 0mm)
  
  gr(size=(500,400))  
end

function plot_scatter(res::Dict, res_prop::Dict; plot_prop=false, size=(800,600))
    L = length(res)
    tmax = length(res[L[1]])
    if plot_prop
      scatter(1:tmax, [res_prop[i] for i in 1:L], xaxis=:log, legendtitle= "level", 
            legend=:outertopright, labels = collect(1:L)', palette = palette(:Blues)[2:(L+1)],
            markerstrokewidth = 0, markersize = 3., size=size)
    else 
      scatter(1:length(res[1]), [res[i] for i in 1:L], xaxis=:log, legendtitle= "level", 
            legend=:outertopright, labels = collect(1:L)', palette = palette(:Reds)[2:(L+1)],
            markerstrokewidth = 0, markersize = 3., size=size)
      global_freq = [sum([res[ℓ][t]*res_prop[ℓ][t] for ℓ in 1:L]) for t in 1:tmax]
      plot!(1:tmax, global_freq, linestyle=:dash, color=:black, width = 1.5, label = "global") 
    end
  end
