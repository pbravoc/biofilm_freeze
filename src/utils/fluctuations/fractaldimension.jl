"""
    _getbins(z::Vector, R::Matrix, S::Vector)
Returns the binning sizes and counts of an array `z`, 
across an array of scaling sizes `S` of a box of dimensions `R`.
"""
function _getbins(z::Vector, R::Matrix, S::Vector)
    rs = R.*S;                                  # Rescale resolution by S
    M=[length(z),maximum(z[isfinite.(z)])];     # Maximums on each dim
    m=[0,minimum(z[isfinite.(z)])];             # Minimums on each dim
    bnum=div.(M-m,rs');                         # Number of bins
    bnum=bnum+iszero.(bnum)                     # making all the bnum=0 bin become 1
    bsiz=(M-m)./bnum                            # adjust bin size
    return bnum',bsiz'
end

"""
    _getcounts(z::Vector, R::Matrix,S::Vector)
Returns the number of unique boxes of the vector `z` after discretization 
of the space across an array of scaling sizes `S` of a box of dimensions `R`.
"""
function getcounts(z::Vector, R::Matrix, S::Vector)
    x = 1:length(z)                 # Discretization on x-dim
    bn, bs = _getbins(z, R, S)      # Get the bin sizes
    counts = zeros(length(S))       # Pre-allocate output
    for i=1:length(S)               # Loop over different sizes
        c=[bn[i,2],1];              # Constants for unique box id.
        C=c[1].*(div.(x,bs[i,1]'))+c[2].*(div.(z,bs[i,2]'))     # Assigning box id.
        C = C'                      # Fix shape orientation
        idx=isfinite.(C);           # Removing NaNs
        counts[i]=length(unique(C[idx]))  # Save unique box count
    end
    return counts
end

"""
    auto_fd(z::Vector;
            z_resolution::Float64 = 1e-3)
Calculates the fractal dimension of a profile, using box-counting. It calculates the 
box-size distribution, and the corresponding counts. Then returns the slope of the  
largest region when fitting in log-log.
"""
function auto_fd(z::Vector; 
                 z_resolution = 1e-3)
    S = 10 .^ Array(0:0.05:log10(length(z)))    # From 1 to the size of z, logspace
    C = getcounts(z, [1 z_resolution], S)      # Get number of unique boxes
    S, C = smoothlog(S, C, 0.3)                 # Smooth the data for better fits
    lr = linear_regions(log10.(S), log10.(C), tol=0.5)      # Find regions in log-log 
    regions_om = lr[1][2:end]-lr[1][1:end-1]    # Find size of large regions 
    idx = sortperm(regions_om, rev=true)[1]     # Index of largest region 
    return round(-lr[2][idx], digits=3)         # Return the slope of the largest region
end

"""
    fixed_fd(z::Vector,
             x_low::Float64,
             x_high::Float64;
             z_resolution::Float64 = 1e-3)
Calculates the fractal dimension of a profile, using box-counting. It calculates the 
box-size distribution, and the corresponding counts. Then returns the slope of the  
region between s_low and s_high in log-log.
"""
function fixed_fd(z::Vector,
                  s_low::Float64=1e1,
                  s_high::Float64=1e3;
                  z_resolution::Float64=1e-3)
    S = 10 .^ Array(0:0.05:log10(length(z)))    # From 1 to the size of z, logspace
    C = getcounts(z, [1 z_resolution], S)      # Get number of unique boxes
    return round(logfit(S, C, s_low, s_high)[1], digits=3)  # Fit in log-log
end