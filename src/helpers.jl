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
    open(path, "a") do io
        for t=1:length(sol.u), ℓ=1:L
            for val in sol.u[t].x[ℓ]
                write(io, "$(t) $(ℓ) $(round(val, digits=14))\n")
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


# processing_sol1(x, n) = sum((collect(0:(n-1)) / n) .* x) / sum(x) 

# function parse_sol(s::String)
    
#     sol = CSV.read(s, DataFrame; header=["timestep", "L", "value"])
#     tmax = maximum(unique(sol.timestep))
#     L = 6
#     inst_level = Dict()
#     inst_level_prop = Dict()
#     lower_limit, upper_limit = 1, 21
#     for t=1:tmax
#         for ℓ=1:L
#             myrange = UnitRange(lower_limit:upper_limit)
#             n = length(sol.value[myrange])
#             x = sol.value[myrange]
#             out = processing_sol1(x, n)
#             out_prop = sum(x)
#         if haskey(inst_level, ℓ)
#             inst_level[ℓ] = [inst_level[ℓ]; out]
#             inst_level_prop[ℓ] = [inst_level_prop[ℓ]; out_prop]
#         else
#             inst_level[ℓ] = out
#             inst_level_prop[ℓ] = out_prop
#         end
  
#         lower_limit += 21
#         upper_limit += 21
  
#         end
#     end
#     return inst_level, inst_level_prop
#   end
  

# function parse_sol(s::ODESolution)
#     L = length(s.u[1].x)
#     tmax = length(s)-1
#     inst_level = Dict()
#     inst_level_prop = Dict()
#     for ℓ=1:L
#       values = []
#       values_prop = []
#       for t=1:tmax
#         n = length(s.u[t].x[ℓ])
#         x = s.u[t].x[ℓ]
#         out = processing_sol1(x,n)
#         push!(values, out)
#         out = sum(x)
#         push!(values_prop, out)
#       end
#       inst_level[ℓ] = values
#       inst_level_prop[ℓ] = values_prop
#     end
#     return inst_level, inst_level_prop
# end

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

# using SQLite

# function plot_phase_diagram(res_db)
    
#     L = length(res)
#     ys = [last(inst_level[i]) for i=1:L]

#     ps = [heatmap(1:length(res[1]), ys[lvl]) for lvl=1:6]
# end