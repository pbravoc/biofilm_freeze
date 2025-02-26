using Agents
using CellListMap
using StaticArrays
using LinearAlgebra 
using Distributions 
using DataFrames 
using JLD2 
using Distances 
using Loess
using Dates 

##
@agent struct Particle(ContinuousAgent{2,Float64})
    r::Float64              # radius
    mass::Float64           # repulsion force constant
    net_force::Float64      # magnitude of the net force, piconewtons 
    θ::Float64              # Direction of the displacement, polar coordinates
    E::Float64              # Elastic modulus
end
Bacteria(; vel, r, mass, net_force, θ, E) = (vel, r, mass, net_force, θ, E)

"""
    initialize_model(seed_height;
                          number_of_particles = 10_000,
                          sides = SVector(50.0, 50.0),
                          dt = 0.001,
                          parallel = true,
                          E = 100.0,
                          r = 0.5,
                          σᵣ = 0.05,
                          max_radius = 4.0)  
                          
Initializes the interacting agent particle model, as well as the underlying space for 
fast force calculations using CellListMap. 

Given the geometry of the system, delimited by `sides`, there's an additional `seed_height`,
that determines the initial binding box of the system, just so it can then relax in a 
physically relevant way.
"""
function initialize_model(seed_height;
                          number_of_particles = 10_000,
                          sides = SVector(50.0, 50.0),
                          dt = 0.001,
                          parallel = true,
                          E = 100.0,
                          r = 0.5,
                          σᵣ = 0.05,
                          max_radius = 4.0)  

    # Start cells randomly in the seeding box
    positions = [SVector{2, Float64}([sides[1], seed_height-0.5].*rand(2) + [0.0, 0.5])
                 for _ in 1:number_of_particles]

    forces = similar(positions)                             # Initialize forces array
    space2d = ContinuousSpace(sides; periodic = false)      # Initialize space

    # Set CellListMap.jl system and link to Agents.jl
    system = ParticleSystem(positions = positions, unitcell = sides, 
                            cutoff = 2*max_radius, output = forces, 
                            output_name = :forces, parallel = parallel)
    properties = (dt = dt, number_of_particles = number_of_particles, system = system)
    model = StandardABM(Particle, space2d;
                        agent_step!, model_step!, 
                        agents_first=false, 
                        properties = properties)

    # Create active agents
    for id in 1:number_of_particles
        new_radius = rand(Normal(r, σᵣ))            # Sample from radius distr
        pos = positions[id]
        new_particle = Bacteria(r = new_radius,                 # Radius [μm]
                                mass = 0.33*π*new_radius^2,     # Weight [pg]
                                net_force = 0.0,                # Force [μN]
                                θ = 0.0,                        # Angle 
                                vel = SVector{2}(0.0, 0.0),     # Velocity
                                E = E)                          # Elastic modulus                 
        add_agent!(pos, Particle, model, new_particle...)       # Add to the model
    end
    return model
end

"""
    calc_forces!(x, y, i, j, d2, forces, model)
Use CellListMap to calculate the forces in the model.
"""
function calc_forces!(x, y, i, j, d2, forces, model)
    pᵢ = model[i]                   # Cell i
    pⱼ = model[j]                   # Cell j
    distance = sqrt(d2)             # Distance between cell centers
    if distance ≤ (pᵢ.r + pⱼ.r)     # If the cells are overlapping
        d = pᵢ.r + pⱼ.r             # Effective diameter
        dr = (y - x)/d              # Normalized distance from i->j (as a vector)
        h = (pᵢ.r + pⱼ.r)^2 - d2    # Proxy for overlap 
        fij = - pᵢ.E * (d^0.5) * (h^1.5) * dr        # Hertzian theory
        forces[i] += fij
        forces[j] -= fij
    end
    return forces
end

"""
    model_step!(model::ABM)
Step all the forces for the particles.
"""
function model_step!(model::ABM)
    # Update the pairwise forces at this step
    map_pairwise!((x, y, i, j, d2, forces) -> calc_forces!(x, y, i, j, d2, forces, model), 
                  model.system)
    return nothing
end

"""
    agent_step!(agent, model::ABM)
Steps for the agent. At the moment this consists of:
    1. Grab forces from CellListMap
    2. Save net forces 
    3. Update positions 
"""
function agent_step!(agent, model::ABM)
    id = agent.id
    dt = abmproperties(model).dt
    f = model.system.forces[id]             # Retrieve the forces on agent id
    agent.net_force = norm(f)               # Save the net force
    agent.θ = dot(f, [0, 1]) / norm(f)      # Get angle
    a = f / agent.mass                      # Get the accelerations
    v = a * dt                              # Overdamped, no velocity 'memory'
    x = agent.pos + v*dt                    # Update position
    x = normalize_position(x, model)        # Normalize positions for boundary
    agent.vel = v                           # Save velocity
    move_agent!(agent, x, model)            # Move agent according to normalized x
    model.system.positions[id] = SVector(agent.pos) # Update positions in the force calc
    return nothing
end

"""
    add_perturbation!(per_radius::Float64, new_position::Tuple, model::ABM)
It adds a new particle of radius `per_radius` at a `new_position`, while correctly 
scaling all the underlying structures that are needed to keep the system running. 
"""
function add_perturbation!(per_radius::Float64, E::Float64, new_position::Tuple, model::ABM)
    new_position = SVector(new_position)
    perturbation = Bacteria(r = per_radius,                 # Radius [μm]
                            mass = 0.33*π*per_radius^2,     # Weight [pg]
                            net_force = 0.0,                # Force [pN]
                            θ = 0.0,                        # Angle 
                            vel = SVector{2}(0.0, 0.0),     # Velocity
                            E = E)          # Elastic modulus           
    add_agent!(new_position, Particle, model, perturbation...)  # Add to the model    
    push!(model.system.positions, new_position)             # Add positions to system 
    resize_output!(model.system, length(model.system.positions))
end

##
""" 
    smooth_topography(h_profile::Tuple;
                     smooth_window::Float64=0.173)
Use a local regression to smooth the topography. We utilize the pixel resolution in the 
Zygo ZeGage Pro interferometer as the smooth_window:
- 5.5x = 1.572 μm/px 
- 10x = 0.865 μm/px 
- 50x = 0.173 μm/px 
"""
function smooth_topography(h_profile::Tuple;
                           smooth_window::Float64=0.173)
    vx = h_profile[1]                                       # Get x values
    vy = h_profile[2]                                       # Get y values 
    model = loess(vx, vy, span=smooth_window/maximum(vx))   # Fit Loess algorithm 
    ny = predict(model, vx)                                 # Use model to get new y values 
    return (vx, ny)
end

"""
    topography(df::DataFrame, 
               padding_distance::Float64; 
               topography_percentile::Float64 = 0.8, 
               dx::Float64 = 1e-2, 
               x_floor::Float64 = 0.0, 
               x_ceil::Float64 = 20.0,
               smooth_window::Float64=0.173)
From the cell positions and sizes defined in `df`, calculate the interface topography 
with a `padding_distance` accounting for Extracellular Polymeric Substances. This function 
works by setting an underlying grid of reference points, and then calculating the distances 
to those points to the closest cell+padding. 
"""
function topography(x_vec::Vector,
                    y_vec::Vector,
                    r_vec::Vector, 
                    padding_distance::Float64; 
                    topography_percentile::Float64 = 0.8, 
                    dx::Float64 = 1e-2, 
                    x_floor::Float64 = 0.0, 
                    x_ceil::Float64 = 50.0,
                    smooth_window::Float64=0.173)
    df = DataFrame()
    df.x = x_vec 
    df.y = y_vec 
    df.r = r_vec
    topo_cutoff = quantile(df.y, topography_percentile) # Get height cutoff
    topo_max = maximum(df.y) + maximum(df.r) + padding_distance 
    topo_cells = filter(x->x.y > topo_cutoff, df)       # Get only cells near the top

    # Generate reference points
    x_list = (x_floor:dx:x_ceil)
    y_list = (topo_cutoff:dx:topo_max)
    my_points = reduce(vcat, [x for x in Iterators.product(x_list, y_list)])
    x_ref, y_ref = [p[1] for p in my_points], [p[2] for p in my_points]

    # Calculate all the distances
    pw = pairwise(PeriodicEuclidean([x_ceil Inf]),  # x boundary conditions at the ceiling
                  [x_ref y_ref]', [topo_cells.x topo_cells.y]', dims=2)

    idx_occupied = falses(size(pw)[1])      # Start boolean index for grid points
    for i=1:size(topo_cells)[1]
        idx_cell = pw[:,i] .< (topo_cells.r[i] + padding_distance)
        idx_occupied = idx_occupied .| idx_cell     # Mark points that are occupied
    end

    reference = DataFrame("x"=>x_ref, "y"=>y_ref, "occupied"=>idx_occupied)
    reference = filter(x->x.occupied == true, reference)    # Select only 'close' references
    gf = groupby(reference, :x)
    gf_result = combine(gf, :y=>maximum)                    # Tallest reference 
    sort!(gf_result, [:x])
    h_smooth = smooth_topography((gf_result.x, gf_result.y_maximum); 
               smooth_window=smooth_window)
    return h_smooth
end

function model_change(model, df_initial; L::Float64 = 50.0)
    df_current = collect_agent_data!(base_df(), model, adata)   # Save current positions
    tf = DataFrame()
    tf.id = df_initial.id
    tf.r = df_initial.r
    tf.x0 = df_initial.x 
    tf.y0 = df_initial.y 
    tf.xf = df_current.x 
    tf.yf = df_current.y 
    tf.mass = df_current.mass
    Δx_tmp = reduce(hcat, [sign.(tf.xf .- tf.x0) .* (tf.xf .- tf.x0),
                                   L .- sign.(tf.xf .- tf.x0) .* (tf.xf .- tf.x0)])
    tf.Δx = [minimum(Δx_tmp[i,:]) for i=1:size(Δx_tmp)[1]]
    tf.Δy = abs.(tf.yf - tf.y0)
    tf.Δr = sqrt.(tf.Δx .^ 2 + tf.Δy .^2 )
    dx_tmp = reduce(hcat, [sign.(tf.x0 .- tf.x0[end]) .* (tf.x0 .- tf.x0[end]),
                           L .- sign.(tf.x0 .- tf.x0[end]) .* (tf.x0 .- tf.x0[end])])
    tf.dx = [minimum(dx_tmp[i,:]) for i=1:size(dx_tmp)[1]]
    tf.dy = tf.y0 .- tf.y0[end]
    tf.dr = sqrt.(tf.dx .^ 2 + tf.dy .^2 )
    tf.height_strat = fld.(tf.y0, 5.0) 
    return tf
end

function perturb_model(file_name, 
                       perturb_strength,
                       perturb_E, 
                       perturb_position; 
                       sim_steps::Int = 1_000)
    model = jldopen(file_name)["model"]                             # Load model 
    add_perturbation!(perturb_strength, perturb_E, perturb_position, model)    # Add perturbation 
    df_i = collect_agent_data!(base_df(), model, adata)             # Save initial positions
    Agents.step!(model, sim_steps)                                  # Step the model
    tf = model_change(model, df_i)                          # Generate comparison dataframe 
    tf.perturb_strengh .= perturb_strength 
    tf.perturb_x .= perturb_position[1]
    tf.perturb_y .= perturb_position[2]
    return tf 
end 

x(agent) = agent.pos[1]     # Get x position
y(agent) = agent.pos[2]     # Get y position
adata = [x, y, :r, :mass, :net_force, :θ]   # Data structure for particle retrieving
base_df() =  DataFrame("time"=>[], "id"=>[], "x"=>[], "y"=>[],  # DataFrame for simulations
                        "r"=>[], "mass"=>[], "net_force"=>[], "θ"=>[])

function optimize_datatypes!(df::DataFrame)
    for col_name in names(df)
        col = df[!, col_name]
        col_type = eltype(col)
        
        # Skip optimization for certain column types
        if col_type <: AbstractString || col_type <: Symbol || col_type <: Date || col_type <: DateTime
            println("Skipping column '$col_name' of type $col_type (already optimized type)")
            continue
        end
        
        # Try to infer the best type
        try
            # Check if all values are missing or nothing
            if all(ismissing, col) || all(isnothing, col)
                println("Column '$col_name' contains only missing/nothing values")
                continue
            end
            
            # For numeric columns
            if col_type <: Number || col_type <: Any
                # Check if all values can be integers
                non_missing = filter(x -> !ismissing(x) && !isnothing(x), col)
                
                if isempty(non_missing)
                    continue
                end
                
                # Check if all values are integers
                all_ints = all(x -> isa(x, Number) && isinteger(Float64(x)), non_missing)
                
                if all_ints
                    # Check the range to determine the most compact integer type
                    min_val = minimum(Float64.(non_missing))
                    max_val = maximum(Float64.(non_missing))
                    
                    if min_val >= 0
                        if max_val <= typemax(UInt8)
                            df[!, col_name] = convert.(UInt8, col)
                            println("Converted column '$col_name' to UInt8")
                        elseif max_val <= typemax(UInt16)
                            df[!, col_name] = convert.(UInt16, col)
                            println("Converted column '$col_name' to UInt16")
                        elseif max_val <= typemax(UInt32)
                            df[!, col_name] = convert.(UInt32, col)
                            println("Converted column '$col_name' to UInt32")
                        else
                            df[!, col_name] = convert.(UInt64, col)
                            println("Converted column '$col_name' to UInt64")
                        end
                    else
                        if min_val >= typemin(Int8) && max_val <= typemax(Int8)
                            df[!, col_name] = convert.(Int8, col)
                            println("Converted column '$col_name' to Int8")
                        elseif min_val >= typemin(Int16) && max_val <= typemax(Int16)
                            df[!, col_name] = convert.(Int16, col)
                            println("Converted column '$col_name' to Int16")
                        elseif min_val >= typemin(Int32) && max_val <= typemax(Int32)
                            df[!, col_name] = convert.(Int32, col)
                            println("Converted column '$col_name' to Int32")
                        else
                            df[!, col_name] = convert.(Int64, col)
                            println("Converted column '$col_name' to Int64")
                        end
                    end
                else
                    # For floating-point numbers, use Float32 to save memory
                    df[!, col_name] = convert.(Float32, col)
                    println("Converted column '$col_name' to Float32")
                end
            end
            
            # For string-like Any columns
            if col_type <: Any
                all_strings = all(x -> ismissing(x) || isnothing(x) || isa(x, AbstractString), col)
                if all_strings
                    df[!, col_name] = convert.(Union{String, Missing}, col)
                    println("Converted column '$col_name' to String")
                end
            end
            
        catch e
            println("Error optimizing column '$col_name': $e")
        end
    end
    return df
end
                        
