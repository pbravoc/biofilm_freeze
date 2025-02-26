using Arrow, DataFrames, DataFramesMeta
using Plots, StatsPlots 
using Statistics
using FractalDimensions 
include("../scripts/FluctAnalysis.jl")
using .FluctAnalysis

function _homeland(profile; hl_buffer = 500)
    half_length = Int((2000 ÷ 0.173) ÷ 2)+hl_buffer   # Homeland length / Resolution / 2
    profile_center = length(profile) ÷ 2        # Center of the profile
    homeland = profile[profile_center-half_length:profile_center+half_length]  
    return homeland 
end

function get_topography(profile;
                        hl_buffer = 500,
                        λ = 500.0,
                        n = 4) 
    # Get homeland range  
    topography = _homeland(profile)
    topography = fillnans(topography)        # Fill NaNs in homeland 
    topography = butterworth_filter(topography, λ; n=n)[hl_buffer:end-hl_buffer] 
    return topography
end

function avg_height(profile)
    homeland = _homeland(profile)
    return mean(homeland)
end

function roughness_calc(topo)
    S_array, w_array = getwloc(topo, 0.173)
    l, w, h = w_data_extraction(S_array, w_array)    # Extract the scalars
    return S_array, w_array, l, w, h
end 

##
df = DataFrame(Arrow.Table("data/vgd_database.arrow"))                  # Load database
df = filter(x->x.replicate in ["A", "B", "C"] && x.time < 48.0 , df)    # Filter data
df.topo = [get_topography(row.profile) for row in eachrow(df)]          # Get topography
df.avg_height = [avg_height(row.profile) for row in eachrow(df)]
##
df = @eachrow df begin                                             # Get ISO Metrics 
    @newcol :Ra::Vector{Float32}
    @newcol :Rz::Vector{Float32}
    @newcol :Rp::Vector{Float32}
    @newcol :Rv::Vector{Float32}
    @newcol :Rq::Vector{Float32}
    @newcol :Rsk::Vector{Float32}
    @newcol :Rku::Vector{Float32}
    :Ra = Ra(:topo)
    :Rz = Rz(:topo)
    :Rp = Rp(:topo)
    :Rv = Rv(:topo)
    :Rq = Rq(:topo)
    :Rsk = Rsk(:topo)
    :Rku = Rku(:topo)
end

df = @eachrow df begin                                          # Power Spectrum 
    @newcol :freq::Vector{Vector{Float32}}
    @newcol :psd::Vector{Vector{Float32}}
    @newcol :psd_exponent::Vector{Float32}
    :freq, :psd, :psd_exponent = get_ν(:topo, 0.173, 500.0)
end

df = @eachrow df begin                                          # Wsat and Hurst
    @newcol :l::Vector{Vector{Float32}}
    @newcol :wl::Vector{Vector{Float32}}
    @newcol :length_cross::Vector{Float32}
    @newcol :w_sat::Vector{Float32}
    @newcol :hurst::Vector{Float32}
    :l, :wl, :length_cross, :w_sat, :hurst = roughness_calc(:topo)
end

##
df = @eachrow df begin 
    @newcol :autocorr::Vector{Vector{Float32}}
    @newcol :autocorr_length::Vector{Float32}
    :autocorr, :autocorr_length = autocorr_length(:topo)
end
##
df = @eachrow df begin                                          # Fractal Dimension
    @newcol :FD::Vector{Float32}
    x = Array(1:length(:topo))
    X = standardize(StateSpaceSet([x :topo]))
    :FD = grassberger_proccacia_dim(X)
end
# Get correlation lengths 

# Save full dataset 

# Save scalar dataset 