include("transforms.jl")
export smoothlog, subpoly, logfit, fillnans

include("linear_regions.jl")
export linear_regions, linear_region

include("timelapsealign.jl")
export timelapsealign!

include("autocorrelation.jl")
export autocorr, autocorr_length

include("background_substract.jl")
export surface_profile, find_edge, sub_background