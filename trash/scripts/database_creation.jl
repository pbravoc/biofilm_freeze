using DataFrames, CSV, Glob 
using Arrow 

"""
    df_vgd_folder(strain_name::String)
    This function returns a dataframe from the vertical growth dynamics paper Dryad `.csv` 
    files. The data for each timelapse is:
        1. Replicate (A, B, C) for each strain 
        2. Time in hours respective to the inoculation 
        3. And the raw height data for the biofilm surface. 
"""
function add_vgd_folder(strain_name::String)
    dir = "data/vertical_growth_dynamics/$strain_name/"     # Directory with strains
    df = DataFrame(CSV.File(dir*"$strain_name.csv"))
    df = filter(x -> x.replicate in ["A", "B", "C"] , df)   # Remove controls and longtime
    N = size(df, 1) / 3
    df.order = Int.(repeat(Array(1:N), 3))
    df = filter(x -> x.order < N, df)
    profile_storage = []
    for replicate in ["A", "B", "C"]
        data = CSV.File(dir*"$(strain_name)_$replicate.csv") |> Tables.matrix;
        profiles = [data[i, :] for i in 1:size(data, 1)]
        push!(profile_storage, profiles)
    end 
    df.profile = reduce(vcat, profile_storage)
    df.strain = repeat([strain_name], size(df)[1])
    select!(df, [:strain, :replicate, :order, :time, :zoom, :profile])
    return df
end

strain_names = ["bgt127", "bh1514", "ea387", "jt305", "sw519", "sw520"]
data = []

for strain in strain_names                      # Loop over different experiments
    strain_dataframe = add_vgd_folder(strain)
    push!(data, strain_dataframe)
end

data = vcat(data...)                            # Aggregate into one file 
##
Arrow.write("data/vgd_database_uncompressed.arrow", data)
##
#CSV.write("data/vgd_database.csv", data)        # Save to a large file 
##
data.length = [length(x) for x in data.profile]
##
plot(data.length, color=:red)