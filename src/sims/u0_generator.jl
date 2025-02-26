include("simulation_functions.jl")

# Initialize model
N_range = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
height_range = [7.0, 14.0, 21.0, 25.0, 30.0, 35.0, 37.0, 40.0]

i = parse(Int, ARGS[1])
j = parse(Int, ARGS[2])

seed_height = height_range[i]
number_of_particles = N_range[i]
folder_loc = "/localdata/Pablo/biofilm_freeze/data/simulations/00_initialcond/"
fname = folder_loc*"u0_"*lpad(string(number_of_particles), 4, "0")*"_"*string(j)*".jld2"
model = initialize_model(seed_height, 
                        number_of_particles = number_of_particles, 
                        sides= SVector(50.0, 100.0))
Agents.step!(model, 10_000_000)              # Very long simulations to ensure equilibrium
save(fname, Dict("model"=>model))
