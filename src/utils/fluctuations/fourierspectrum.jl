using FFTW
using DSP 

"""
    powerspectrum(y::Vector, x_space::Float64=1.0)
Returns the power spectrum of the vector `v`, obtaining the fourier
transform and then squaring it.
"""
function powerspectrum(y::Vector, x_space::Float64=1.0)
    y = fillnans(y)                                # Remove NaNs
    fourier_transform = fft(y)                      # Actual fourier transform
    freqs = fftfreq(length(y), 1/x_space)           # Frequencies
    idx = freqs .> 0                                # Get positive frequencies
    freqs = freqs[idx]
    power = (abs.(fourier_transform[idx])) .^ 2     # 
    return freqs, power
end

""" 
    decompose_profile(measured_profile::Vector; 
                      rpx::Float64=0.173,
                      λs::Float64=0.5, 
                      λc::Float64=300.0, 
                      λf::Float64=1000)
Decomposes the `measured_profile` according to the 3 characteristic 
wavelengths `λs`, `λc`, and `λf`
Returns the primary, waviness, and roughness profiles. 
REMOVING IN THE FUTURE. Use butterworth_filter and functions from DSP.jl instead.
"""
function decompose_profile(measured_profile; 
                           rpx::Float64=0.173,      
                           λs::Float64=0.5, 
                           λc::Float64=300.0, 
                           λf::Float64=1000.0)
    measured_profile = fillnans(measured_profile) 
    x = Array(1:length(measured_profile))*rpx
    fourier_transform = fft(measured_profile)
    freqs = fftfreq(length(x), 1/(x[2]-x[1]))
    wavelengths = 1 ./ abs.(freqs)
    primary_idx = wavelengths .> λs
    roughness_idx = (wavelengths .> λs) .* (wavelengths .< λc)
    waviness_idx = (wavelengths .> λc) .* (wavelengths .< λf)
    fourier_primary = copy(fourier_transform)
    fourier_roughness = copy(fourier_transform)
    fourier_waviness = copy(fourier_transform)
    fourier_primary[primary_idx] .= 0 
    fourier_roughness[roughness_idx] .= 0
    fourier_waviness[waviness_idx] .= 0 
    p_profile = measured_profile - real(ifft(fourier_primary))
    r_profile = measured_profile - real(ifft(fourier_roughness))
    w_profile = measured_profile - real(ifft(fourier_waviness))
    return p_profile, w_profile, r_profile
end

""" 
    get_ν(z::Vector; 
          k::Integer==3, 
          rpx::Float64=0.173,
          smooth_window = 0.1,
          tol=0.3)
Gets the ν scaling exponent for the linear regions in the 
log-log power spectrum of the profile `z`. After obtaining the 
spectrum, it is smoothed through LOESS, and then re-sampled such 
that the frequencies are homogeneous in log-space. This ensures 
a proper weighting of the scaling law. `k` is an optional argument
that returns how many of the largest linear regions should be taken 
in account for the calculation.
"""
function get_ν(z::Vector,
               λ_low, 
               λ_high; 
               k = 3, 
               rpx = 0.173, 
               smooth_window = 0.1, 
               tol=0.3)
    f, p = powerspectrum(z, rpx)
    fs, ps = smoothlog(f, p, smooth_window)
    idx = (fs .> λ_low) .* (fs .< λ_high)
    fs = fs[idx]
    ps = ps[idx]
    lr = linear_regions(log10.(fs), log10.(ps), tol=tol)
    ν = _weigh_regions(lr, k)
    return f, p, -ν
end


""" 
    butterworth_filter(z::Vector, λ0::Float64; 
                       λs=0.173, n=4)
Filter the raw profile using a Butterworthfitler of order `n`, cutting off 
wavelengths larger than λ0. Returns a vector of the same size as `z`.
"""
function butterworth_filter(z::Vector, λ0::Float64; 
                            λs=0.173, n=4)
    responsetype = Highpass(1/λ0, fs=1/λs)
    designmethod = Butterworth(n)
    but_filter = digitalfilter(responsetype, designmethod)
    return filtfilt(but_filter, z)
end
