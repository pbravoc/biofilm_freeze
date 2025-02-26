include("localwidth.jl")
export getwloc, w_data_extraction, df_hurst!

include("profile_elements.jl")
export find_profile_elements, element_counter

include("surfaceroughness.jl")
export Ra, Rp, Rv, Rz, Rc, Rq, Rsk, Rku, RSm

include("fourierspectrum.jl")
export powerspectrum, decompose_profile, get_ν, butterworth_filter
 
include("fractaldimension.jl")
export getcounts, auto_fd, fixed_fd