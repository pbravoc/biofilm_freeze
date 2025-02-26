using Polynomials 
using Loess 
using Interpolations

##
"""
    Returns the average (in logspace) between number `a` and `b`.
"""
function _logspace_average(a::Float64, b::Float64)
    a_log = log10(a)
    b_log = log10(b)
    return 10^(0.5*a_log+0.5*b_log)
end

"""
    smoothlog(x::Vector, 
              y::Vector; 
              smooth_window::Float64=0.1,
              x_resolution::Float64=5e-2)
Returns a smoothed out, equally spaced(`x_resolution`) in log-space.
Input is the vectors `x` and `y`. Given the nature of the transform 
the size of the input is NOT the same as the output.
"""
function smoothlog(x, y, smooth_window =0.1, x_resolution =5e-2)
  model = loess(log10.(x), log10.(y), span=smooth_window)       # Loess for smoothing
  x_homogeneous = Array(log10(x[1]):x_resolution:log10(x[end])) # Sample points
  y_homogeneous = predict(model, x_homogeneous)                 # Measure from sample
  return 10 .^ x_homogeneous, 10 .^ y_homogeneous               # Return in linear space
end 

"""
    Handy fitting tool, taking a collection of linear regions, 
and then weighting the slope values depending on the size range
"""
function _weigh_regions(lr::Tuple, k::Integer=3)
    region_sizes = lr[1][2:end] - lr[1][1:end-1]
    idx = sortperm(region_sizes, rev=true)          # Index of the largest regions
    size_sorted = region_sizes[idx][1:k]            # Sort all regions    
    slope_sorted = lr[2][idx][1:k]                  # As well as slopes
    weighted_slope = sum(size_sorted .* slope_sorted) / sum(size_sorted)
    return weighted_slope
end

"""
    fillnans(y::Vector)
Interpolation to remove `NaNs` in a profile before any processing. 
It simply uses a linear interpolation between the surrounding points WITH data. 
Thus, it assumes the data has no **big jumps**, or any interesting features in 
said region. 

*Remember: good measurements are better than any post-processing you can do.*
"""
function fillnans(y::Vector)
    x = Array(1:length(y))          # Get the x values
    idx = (isnan.(y))               # Values that are NaNs
    if sum(idx) < length(y)           # If there are manageable NaNs 
        real_valued = findall(x->x==false, idx) # Get the real values
        start_point = real_valued[1]
        end_point = real_valued[end]
        itp = linear_interpolation(x[.! idx],   # No extrapolation
                                   y[.!idx])
        return itp(x[start_point:end_point])          
    else 
        return y
    end
end

""" 
    subpoly(y::Vector, deg::Int)
Returns the array `y` after substracting the best-fit polynomial of degree `deg`.
"""
function subpoly(y::Vector, deg::Int)
    x = Array(range(1, length(y), step=1))
    valids = .~isnan.(y)
    flat = Polynomials.fit(x[valids], y[valids], deg) # degree = 2
    myfit = flat.(x)
    sub = y-myfit
    return sub
end

"""
    shifted_log(y::Vector, shift::Float64)
Shifts the vector in a quantity `shift` in log-space. This function 
is solely for the purpose of visualization multiple spectrums, akin 
to ridgeplots.
"""
function _shifted_log(y::Vector, shift::Float64)
    ys = log10.(y)      # Get the log space
    ys .+= shift        # Shift
    return 10 .^ (ys)   # Return to linear
end 

"""
    logfit(x::Vector,
            y::Vector,
            x_low::Float64,
            x_high::Flota64)
Fits a set of two arrays `x` and `y`, in the range `[x_low, x_high]`. 
"""
function logfit(x::Vector, 
                y::Vector,
                x_low::Float64,
                x_high::Float64)
    xs, ys = smoothlog(x, y)                        # Smooth for noise removal
    idx = (xs .> x_low) .* (xs .< x_high)           # Indices in the range
    f = fit(log10.(xs[idx]), log10.(ys[idx]), 1)    # Fit using Polynomials
    return -f.coeffs[2], f.coeffs[1]               # Return only the slope
end