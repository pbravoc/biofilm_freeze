using DataFrames 
using NaNMath 

_rms(y1, y2) = sqrt(NaNMath.mean((y1 .- y2).^2 ))

function _aligner(row_1, row_2, lims, n, dx=Array(-600:3:600))
    v1 = subpoly(row_1.profile[lims[1]:lims[2]], n)
    deltas = zeros(length(dx))
    for i in 1:length(dx)
        x = dx[i]
        v2 = subpoly(row_2.profile[lims[1]+x:lims[2]+x], n)
        deltas[i] = _rms(v1, v2)
    end
    idx = findmin(deltas)[2]
    return dx[idx]
end

"""
    timelapsealign!(df::DataFrame)
Aligns, using the rms between profiles at h(t) and  
"""
function timelapsealign!(df::DataFrame)
      lims = [df.hL[1], df.hR[2]]
      Dx = zeros(size(df)[1])
      for i=2:size(df)[1]
          Dx[i] = _aligner(df[i-1,:], df[i,:], lims, 2) 
      end
      dx = Int.(cumsum(Dx))
      df.hL = df.hL + dx 
      df.hR = df.hR + dx
    end