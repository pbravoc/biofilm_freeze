# Code for linear_regions is taken directly from ChaosTools.jl
# To avoid dependency issues we recycle it.

#export linear_region, linear_regions, estimate_boxsizes, linreg
#####################################################################################
# Functions and methods to deduce linear scaling regions
#####################################################################################
using Statistics
using Statistics: covm, varm
# The following function comes from a version in StatsBase that is now deleted
# StatsBase is copyrighted under the MIT License with
# Copyright (c) 2012-2016: Dahua Lin, Simon Byrne, Andreas Noack, Douglas Bates,
# John Myles White, Simon Kornblith, and other contributors.
"""
    linreg(x, y) -> a, b
Perform a linear regression to find the best coefficients so that the curve:
`z = a + b*x` has the least squared error with `y`.
"""
function linreg(x::AbstractVector, y::AbstractVector)
    # Least squares given
    # Y = a + b*X
    # where
    # b = cov(X, Y)/var(X)
    # a = mean(Y) - b*mean(X)
    if size(x) != size(y)
        throw(DimensionMismatch("x has size $(size(x)) and y has size $(size(y)), " *
            "but these must be the same size"))
    end
    mx = Statistics.mean(x)
    my = Statistics.mean(y)
    # don't need to worry about the scaling (n vs n - 1)
    # since they cancel in the ratio
    b = covm(x, mx, y, my)/varm(x, mx)
    a = my - b*mx
    return a, b
end

slope(x, y) = linreg(x, y)[2]


"""
    linear_regions(x, y; dxi::Int = 1, tol = 0.25) -> (lrs, tangents)
Identify regions where the curve `y(x)` is linear, by scanning the
`x`-axis every `dxi` indices sequentially
(e.g. at `x[1] to x[5], x[5] to x[10], x[10] to x[15]` and so on if `dxi=5`).
If the slope (calculated via linear regression) of a region of width `dxi` is
approximatelly equal to that of the previous region,
within tolerance `tol`,
then these two regions belong to the same linear region.
Return the indices of `x` that correspond to linear regions, `lrs`,
and the _correct_ `tangents` at each region
(obtained via a second linear regression at each accumulated region).
"""
function linear_regions(
        x::AbstractVector, y::AbstractVector;
        method = :sequential, dxi::Int = method == :overlap ? 3 : 1, tol = 0.25,
    )
    @assert length(x) == length(y)
    return if method == :overlap
        linear_regions_overlap(x, y, dxi, tol)
    elseif method == :sequential
        linear_regions_sequential(x, y, dxi, tol)
    end
end

function linear_regions_sequential(x, y, dxi, tol)
    maxit = length(x) ÷ dxi

    tangents = Float64[slope(view(x, 1:max(dxi, 2)), view(y, 1:max(dxi, 2)))]

    prevtang = tangents[1]
    lrs = Int[1] #start of first linear region is always 1
    lastk = 1

    # Start loop over all partitions of `x` into `dxi` intervals:
    for k in 1:maxit-1
        tang = slope(view(x, k*dxi:(k+1)*dxi), view(y, k*dxi:(k+1)*dxi))
        if isapprox(tang, prevtang, rtol=tol, atol = 0)
            # Tanget is similar with initial previous one (based on tolerance)
            continue
        else
            # Tangent is not similar.
            # Push new tangent for a new linear region
            push!(tangents, tang)

            # Set the START of a new linear region
            # which is also the END of the previous linear region
            push!(lrs, k*dxi)
            lastk = k
        end

        # Set new previous tangent (only if it was not the same as current)
        prevtang = tang
    end
    push!(lrs, length(x))
    # create new tangents that do have linear regression weighted
    tangents = Float64[]
    for i in 1:length(lrs)-1
        push!(tangents, linreg(view(x, lrs[i]:lrs[i+1]), view(y ,lrs[i]:lrs[i+1]))[2])
    end
    return lrs, tangents
end

"""
    linear_region(x, y; kwargs...) -> ((ind1, ind2), slope)
Call [`linear_regions`](@ref) and identify and return the largest linear region
and its slope. The region starts and stops at `x[ind1:ind2]`.
The keywords `dxi, tol` are propagated as-is to [`linear_regions`](@ref).
The keyword `ignore_saturation = true` ignores saturation that (sometimes) happens
at the start and end of the curve `y(x)`, where the curve flattens.
The keyword `sat = 0.01` decides what saturation is (while `abs(y[i]-y[i+1])<sat` we 
are in a saturation regime).
The keyword `warning = true` prints a warning if the linear region is less than 1/3
of the available x-axis.
"""
function linear_region(x::AbstractVector, y::AbstractVector;
    dxi::Int = 1, tol::Real = 0.2, ignore_saturation = true, warning = true, sat = 0.01)

    isat = 0
    if ignore_saturation
        j = findfirst(i -> abs(y[i] - y[i-1]) > sat, length(y):-1:2)
        if !isnothing(j)
            i = (length(y):-1:2)[j]
            x, y = x[1:i], y[1:i]
        end
        k = findfirst(i -> abs(y[i+1] - y[i]) > sat, 1:length(y)-1)
        if !isnothing(k)
            x, y = x[k:end], y[k:end]
            isat = k-1
        end
    end

    lrs, tangents = linear_regions(x,y; dxi, tol)
    # Find biggest linear region:
    j = findmax(diff(lrs))[2]
    if lrs[j+1] - lrs[j] ≤ length(x)÷3 && warning
        @warn "Found linear region spans less than a 3rd of the available x-axis "*
              "and might imply inaccurate slope or insufficient data. "*
              "Recommended: plot `x` vs `y`."
    end
    return (lrs[j] + isat, lrs[j+1] + isat), tangents[j]
end


