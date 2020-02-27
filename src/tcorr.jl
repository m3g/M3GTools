#
# This function computes the time correlation function of two arrays
#

tcorr(x :: Vector) = tcorr(x,x)

function tcorr(x_in :: Vector, y_in :: Vector)

  x = copy(x_in)
  y = copy(y_in)

  n = size(x)[1]
  if size(y)[1] != n
    error(" In tcorr both vectors are expected to have the same dimensions. ")
  end

  # Centralizing at the average 

  xav = 0.
  yav = 0.
  @inbounds for i in 1:n
    xav = xav + x[i]
    yav = yav + y[i]
  end
  xav = xav / n
  yav = yav / n
  @inbounds for i in 1:n
    x[i] = x[i] - xav
    y[i] = y[i] - xav
  end
  
  # Computing the time-dependent correlation function

  xsq = Vector{Float64}(undef,n)
  ysq = Vector{Float64}(undef,n)
  @inbounds for i in 1:n
    xsq[i] = x[i]^2
    ysq[i] = y[i]^2
  end

  tcorr = zeros(Float64,n)

  @inbounds for dt in 0:n-1
    xnorm = 0.
    ynorm = 0.
    for i in 1:n-dt
      tcorr[dt+1] = tcorr[dt+1] + x[i]*y[i+dt]
      xnorm = xnorm + xsq[i]
      ynorm = ynorm + ysq[i+dt]
    end
    tcorr[dt+1] = tcorr[dt+1] / sqrt(xnorm*ynorm)
  end

  t = collect(1:n)

  return t, tcorr

end
