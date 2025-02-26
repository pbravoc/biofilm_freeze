#= Surface roughness metrics for 1-dimensional profiles
https://www.keyence.com/ss/products/microscope/roughness/line/parameters.jsp
Metrics from the ISO 4287:1997 standard. 
=#
using NaNMath

_moment(z::Vector, i::Integer) = NaNMath.mean((z .- NaNMath.mean(z)).^i)

Ra(z::Vector) = NaNMath.mean(abs.(z .- NaNMath.mean(z)))

Rz(z::Vector) = NaNMath.maximum(z) - NaNMath.minimum(z)

Rp(z::Vector) = NaNMath.maximum(z) - NaNMath.mean(z)
 
Rv(z::Vector) = abs.(NaNMath.minimum(z) - NaNMath.mean(z))
    
Rq(z::Vector) = sqrt(_moment(z, 2))

Rsk(z::Vector) = _moment(z, 3) / (_moment(z, 2) .^ (1.5))

Rku(z::Vector) = _moment(z, 4) / (_moment(z, 2) .^ (2))

Rc(z::Vector) = mean(Rz.([z[bounds[1]:bounds[2]] for bounds in find_profile_elements(z)]))

RSm(z::Vector) = mean([bounds[2]-bounds[1] for bounds in find_profile_elements(z)])
##