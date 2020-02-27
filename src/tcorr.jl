#
# This function computes the time auto-correlation function of a vector
#

function tcorr(x_in :: Vector{Float64}) 

  x = copy(x_in)

  n = size(x)[1]

  # Centralizing at the average 

  xav = 0.
  @inbounds for i in 1:n
    xav = xav + x[i]
  end
  xav = xav / n
  @inbounds for i in 1:n
    x[i] = x[i] - xav
  end
  
  # Computing the time-dependent correlation function

  tcorr = zeros(Float64,n)

  @inbounds for dt in 0:n-1
    xnorm1 = 0.
    xnorm2 = 0.
    for i in 1:n-dt 
      tcorr[dt+1] = tcorr[dt+1] + x[i]*x[i+dt]
      xnorm1 = xnorm1 + x[i]^2
      xnorm2 = xnorm2 + x[i+dt]^2
    end
    tcorr[dt+1] = tcorr[dt+1] / sqrt(xnorm1*xnorm2)
  end

  for i in 1:n
    x[i] = x[i] + xav
  end

  lags = collect(0:n-1)
  return lags, tcorr

end

#
# This function computes the time correlation function of two vectors
#

function tcorr(x_in :: Vector{Float64}, y_in :: Vector{Float64})

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
  
  tcorr = zeros(Float64,n)

  @inbounds for dt in 0:n-1
    xnorm = 0.
    ynorm = 0.
    for i in 1:n-dt
      tcorr[dt+1] = tcorr[dt+1] + x[i]*y[i+dt]
      xnorm = xnorm + x[i]^2
      ynorm = ynorm + y[i+dt]^2
    end
    tcorr[dt+1] = tcorr[dt+1] / sqrt(xnorm*ynorm)
  end

  lags = collect(0:n-1)
  return lags, tcorr

end
