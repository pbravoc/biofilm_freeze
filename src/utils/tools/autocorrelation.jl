using DSP 

""" 
    autocorr(z::Vector)
Returns the normalized autocorrelation of the vector `z`, since it is symmetric by 
definition, we return only the right side of the function.
"""
function autocorr(z::Vector)
    cross_corr = xcorr(z, z; padmode = :none)       # Centered correlation
    cross_corr = cross_corr[end-length(z)+1:end]    # Get only one side 
    return cross_corr ./ cross_corr[1]              # Return normalized values 
end
 
"""
    autocorr_length(z::Vector;
                    res::Float64 = 0.173,
                    xcorr_value::Float64=0.0)
Returns the autocorrelation length for a vector `z`, which is uniform with a separation 
`res`. The criterion value is `xcorr_value`, by default we set it to 0. Depending on the 
community, it varies, so a clear specification is needed.
"""
function autocorr_length(z::Vector; 
                         res::Float64 = 0.173,
                         xcorr_value::Float64 = 0.0)
    cross = autocorr(z)
    idx = findfirst(x-> x < xcorr_value, cross)
    if idx === nothing 
        return cross, NaN 
    else
        return cross, round(idx*res, digits=3) 
    end
end
