using CircularArrays 
using NaNMath
using DataFrames 

""" 
    _mirror(z::Vector, s::Integer)
Returns a shifted vector `z` by an integer number of spaces `s`. 
NaNs are placed when no shift is allowed
"""
function _mirror(z::Vector, spaces::Integer)
    mirrored_vector = circshift(z, spaces)
    if spaces>=0
        mirrored_vector[1:spaces] .= NaN
    else 
        mirrored_vector[end+spaces+1:end] .= NaN
    end
    return mirrored_vector
end

"""
    _update_moments!(z::Vector, 
                     r::Vector, 
                     counts::Vector, 
                     moment_1::Vector, 
                     moment_2::Vector)
Returns the updated moments for a range of sizes `s`, in a vector `z`. 
The first (`moment_1`), and second (`moment_2`) moments are updated, as well as 
the total number of valid counts `c` at said location.
"""
function _update_moments!(z::Vector, 
                          r::Integer, 
                          counts::Vector, 
                          moment_1::Vector, 
                          moment_2::Vector)
    for current_shift in (r, -r)
        mirrored_vector = _mirror(z, current_shift)
        idx = isnan.(mirrored_vector)                  # values that are NaNs
        counts .+= .! idx                   # Only count values with data
        mirrored_vector[idx] .= 0           # NaNs->0 to not remove data
        moment_1 .+= mirrored_vector        # Update first moment
        moment_2 .+= mirrored_vector .^ 2   # Update second moment
    end
    return counts, moment_1, moment_2
end

"""
    getwloc(z::Vector, 
            rpx::Float;
            rx::Float64=0.4)
Returns the local width of a vector `z`, across an increasing range of sampling size.
Sampling size goes from 1 to the size of `z`. Resolution is given by the optional argument 
`rx` in log2 space.
As the standard deviation is given by:
```math
\\sigma = \\sqrt{E[X^2]-E[X]^2}
```
We keep a cache of the first and second moments as we sample from larger regions.
"""
function getwloc(z::Vector,
                 rpx::Float64; 
                 rx::Float64=0.4)
    S = Int.(floor.(2 .^ (0:rx:log2(length(z)))))   # Points as integer `distances`
    counts = zeros(length(z))                       # Initialize counts and moments
    moment_1 = zeros(length(z))
    moment_2 = zeros(length(z))
    wloc = zeros(length(S)-1)                       # Initialize local widths
    for i in 1:(length(S)-1)                        # Loop over distances
        counts, moment_1, moment_2 = _update_moments!(z,        # Update info
                                                      S[i+1], 
                                                      counts, 
                                                      moment_1, 
                                                      moment_2)
        norm_moment_1 = moment_1 ./ counts          # Normalize first moment
        norm_moment_2 = moment_2 ./ counts          # Normalize second moment
        variances = norm_moment_2 .- (norm_moment_1 .^ 2)   # Calculate variances
        wloc[i] = NaNMath.mean(sqrt.(abs.(variances)))  # Update local width 
    end
    return S[2:end] .* rpx, wloc                    # Return sizes and local widths
end

"""
    wloc_data(distances::Vector,
              local_widths::Vector)
From a characteristic local width profile, extract the 3 key metrics:
1. l_sat: saturation length for the power law.
2. w_sat: saturation width.
3. hurst: scaling exponent in the scaling region.
"""
function w_data_extraction(distances, local_widths)
    spec_x, spec_y = smoothlog(distances, local_widths)         # Smooth data 
    lr = linear_regions(log10.(spec_x), log10.(spec_y), tol=0.4)    # Find linear in log-log
    regions_om = lr[1][2:end]-lr[1][1:end-1]                    # Region size
    idx = sortperm(regions_om, rev=true)[1]                     # Find largest region 
    l_sat = spec_x[sum(regions_om[1:idx])]                      # Saturation length 
    w_sat = spec_y[end]                                         # Saturation width 
    hurst = lr[2][idx]                                          # Hurst exponent
    return l_sat, w_sat, hurst                                  # return all
end

"""
    df_hurst!(df::DataFrame)
Calculate the scaling law for the saturation width `wloc`, and then extract the associated 
scalar metrics `l_sat`, `w_sat`, `hurst`. Adds the quantities to 
"""
function df_hurst!(df::DataFrame)
    df.rpx = [row.zoom .== 50 ? 0.173 : 0.865 for row in eachrow(df)]
    locs, wlocs = [], [] 
    for row in eachrow(df)
        l, w = getwloc(row.homeland, row.rpx; rx=0.3)   # Get scaling arrays
        append!(locs, [l])
        append!(wlocs, [w])
    end
    df.loc = locs
    df.wloc = wlocs
    l_sat = []
    w_sat = []
    hurst = []
    for row in eachrow(df)
        l,w,h = w_data_extraction(row.loc, row.wloc)    # Extract the scalars
        append!(l_sat, [l])
        append!(w_sat, [w])
        append!(hurst, [h])
    end
    df.l_sat = l_sat 
    df.w_sat = w_sat
    df.hurst = hurst
    return df                   # Filter out loc and wloc?
end