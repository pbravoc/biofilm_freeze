using NaNStatistics 
using Loess  
using Peaks 
using ImageFiltering

function surface_profile(data; kernel_size=5)
    patch_size = (kernel_size, kernel_size)
    img = mapwindow(nanmedian, data.surface, patch_size)
	y = nanmean(img, dims=1)' 	# Get average data points 
	z = Vector{Float64}(undef, length(y))
	for i=1:length(y)
		z[i] = y[i]
	end
	idx = .! isnan.(z) 				# Find index with data 
	x = Array(1:length(z))
	loess_model = loess(x[idx], y[idx]; span=0.025)
	z = predict(loess_model, x)
	return subpoly(z, 1)./1000
end

function find_edge(prof)
    x = Array(1:length(prof))
    g2 = (prof[1:end-2] + prof[3:end] - 2*prof[2:end-1])
    smooth_g2 = predict(loess(x[2:end-1], g2; span=0.05), x[2:end-1])
    pks, vals = findmaxima(smooth_g2, 50)
    left_side = floor(Int64, length(smooth_g2)/3)
    right_side = length(smooth_g2)-left_side
    left_idx = argmax(smooth_g2[pks[pks .< left_side]])
    right_idx = argmax(smooth_g2[pks[pks .> right_side]])
    left_border = pks[pks .< left_side][left_idx]
    right_border = pks[pks .> right_side][right_idx]
    return left_border, right_border 
end

function sub_background(y::Vector, le, re, deg::Int)
    x = Array(range(1, length(y), step=1))
    idx = (x .< le - 200) .| (x .> re + 200)
    flat = Polynomials.fit(x[idx], y[idx], deg)
    return y .- flat.(x)
end