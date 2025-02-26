#=An element of a  line profile corresponds to a peak-valley combo.
The requisites to discern between an element and noise are:
    1. Height > 10% maximum Height
    2. Length > 1% total length
=#
using Statistics

function _changeinbool(b::BitVector)
    new_vector = Int.(b[1:end-1]) .+ Int.(b[2:end])
    change_points = findall(x-> x==1, new_vector)
    in_change = change_points[.!b[change_points]]
    out_change = change_points[b[change_points]]
    return in_change, out_change
end

"""
    find_profile_elements(z::Vector)
Returns a list of bounds for the profile elements in a given vector `z`.
A profile element consist of a peak and a valley, which:
    1. rz(peak) >= 0.1 * rz(z)
    2. length(peak) >= 0.01 * length(z) TODO
"""
function find_profile_elements(z::Vector)
    z_min = 0.1*Rz(z)
    peaks = copy(z)
    peak_idx = peaks .> z_min
    valleys = copy(z)
    valley_idx = valleys .< -z_min
    in_valley, out_valley = _changeinbool(valley_idx)
    in_peak, out_peak = _changeinbool(peak_idx)
    transition_points = []
    for i=1:length(out_valley)-1
        left_bound = out_valley[i]
        possible_bounds = findall(x -> left_bound < x, in_peak)
        if !isempty(possible_bounds)
            right_bound = in_peak[possible_bounds[1]]
            if right_bound < out_valley[i+1]
                transition_region = z[left_bound:right_bound]
                zero_idx = findall(x-> abs(x) < 1e-2, transition_region)
                transition_point = left_bound + Int(floor(median(zero_idx)))
                append!(transition_points, transition_point)
            end
        end
    end
    if in_peak[1] < transition_points[1]
        pushfirst!(transition_points, 1)
    else 
        transition_points[1] = 1
    end
    push!(transition_points, length(z))
    element_bounds = [[transition_points[i], transition_points[i+1]] for i=1:length(transition_points)-1]
    return element_bounds
end

"""
    element_counter(z::Vector)
Returns the total amount of elements inside a profile
"""
function element_counter(z::Vector)
    my_elements = find_profile_elements(z)
    return length(my_elements)
end
