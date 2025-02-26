using Plots, StatsPlots
using Glob
include("simulation_functions.jl")

# Initialize model
files = glob("*.jld2", "/localdata/Pablo/biofilm_freeze/data/simulations/00_initialcond/")
model = load(files[10])["model"]
data = collect_agent_data!(base_df(), model, adata)
@df data scatter(:x, :y, zcolor=:net_force, markersize=:mass*20, label=false, title=round(mean(data.net_force), digits=3), clim=(1e-4, 0.2))
plot!(aspect_ratio=:equal, ylim=(0, 60), xlim=(0, 50), grid=false, size=(1200, 820))
#add_perturbation!(0.7, 100.0,  (26.0, 1.5), model)    # Add perturbation 
##
@gif for i=1:100
    Agents.step!(model, 500)              
    data = collect_agent_data!(base_df(), model, adata)
    @df data scatter(:x, :y, zcolor=:net_force, markersize=:mass*20, label=false, title=round(mean(data.net_force), digits=3), clim=(1e-4, 4))
    plot!(aspect_ratio=:equal, ylim=(0, 60), xlim=(0, 50), grid=false, size=(1200, 820))
end
##

#=
N_particles = []
replicate = []
max_height = []
net_force = []
for file in files
    model = load(file)["model"]
    data = collect_agent_data!(base_df(), model, adata)
    data = filter(row->row.y > 3.0, data)
    push!(N_particles, size(data)[1])
    push!(replicate, parse(Int, file[68:71]))
    push!(max_height, maximum(data.y))
    push!(net_force, mean(data.net_force))
end
df = DataFrame("N_particles"=>N_particles, "replicate"=>replicate, "max_height"=>max_height, "net_force"=>net_force)
##
@df df scatter!(:N_particles, :net_force, label="Filter")

i = 35
data = collect_agent_data!(base_df(), load(files[i])["model"], adata)
@df data scatter(:x, :y, zcolor=:net_force, markersize=:mass*20, label=false, title=round(mean(data.net_force), digits=3), clim=(1e-4, 0.2))
plot!(aspect_ratio=:equal, ylim=(0, 60), xlim=(0, 50), grid=false, size=(1200, 820))

=#