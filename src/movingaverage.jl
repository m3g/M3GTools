#
# Function that computes a moving average to smooth plots
#

function movingaverage( x :: Vector, n :: Int )
  
  if ! isodd(n)
    n = n + 1
  end
  delta = round(Int64,(n-1)/2)

  nx = length(x)
  y = similar(x) 

  for i in 1:nx
    y[i] = 0.
    jmin = max(i-delta,1)
    jmax = min(i+delta,nx)
    for j in jmin:jmax
      y[i] = y[i] + x[j]
    end
    y[i] = y[i] / (jmax-jmin+1)

  end

  return y

end

