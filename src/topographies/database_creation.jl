"""
This file aggregates the Vertical Growth Dynamics from Dryad 
(https://datadryad.org/stash/dataset/doi:10.5061/dryad.pg4f4qrsw)
And converts them into a single dataframe that will be used for analysis.
"""

using DataFrames, CSV, Glob 

"""
    df_vgd_folder(strain_name::String)
This function returns a dataframe from the vertical growth dynamics paper Dryad `.csv` 
files. The data for each timelapse is:
    1. Replicate (A, B, C) for each strain 
    2. Time in hours respective to the inoculation 
    3. And the raw height data for the biofilm surface. 
"""
function add_vgd_folder(strain_name::String)
    dir = "data/vertical_growth_dynamics/$strain_name/"
    df = DataFrame(CSV.File(dir*"$strain_name.csv"))
    df = filter(x -> x.replicate in ["A", "B", "C"], df)
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
    select!(df, [:strain, :replicate, :order, :time, :profile])
    return df
end

strain_names = ["bgt127", "bh1514", "ea387", "jt305", "sw519", "sw520"]
data = []

for strain in strain_names
    strain_dataframe = add_vgd_folder(strain)
    push!(data, strain_dataframe)
end

Arrow.write("data/vgd_database.arrow", vcat(data...))
