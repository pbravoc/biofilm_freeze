using Glob, CSV
include("simulation_functions.jl")

i = parse(Int, ARGS[1])
j = parse(Int, ARGS[2])

files = glob("*.jld2", "/localdata/Pablo/biofilm_freeze/data/simulations/00_initialcond/")[6:10]
perturb_folder = "/localdata/Pablo/biofilm_freeze/data/simulations/01_perturbations/"
u0_data = files[i][68:73]
data_file = "perturb_"*u0_data*"_"*lpad(string(j), 3, "0")
output_name = data_file*".csv"
perturb_location = Tuple(rand(2).*(50.0, 17.0))
data = perturb_model(files[i], 0.5, 100.0, perturb_location; sim_steps=1_000_000)
optimize_datatypes!(data)
CSV.write(perturb_folder*output_name, data)