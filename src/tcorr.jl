#
# This function computes the time auto-correlation function of a vector
#

using LinearAlgebra

function tcorr(x_in :: Vector{Float64}) 

  x = copy(x_in)
  n = size(x)[1]

  # Centralizing at the average 

  xav = sum(x)/n
  @. x = x - xav
  
  # Computing the time-dependent correlation function

  tcorr = zeros(Float64,n)
  xnorm1 = 0.
  xnorm2 = 0.
  @inbounds for dt in n-1:-1:0
    xnorm1 = xnorm1 + x[n-dt]^2
    xnorm2 = xnorm2 + x[dt+1]^2
    if xnorm1 != 0 || xnorm2 != 0
      tcorr[dt+1] = LinearAlgebra.dot(@view(x[1:n-dt]),@view(x[dt+1:n]))
      tcorr[dt+1] = 0.
    else
      tcorr[dt+1] = tcorr[dt+1] / sqrt(xnorm1*xnorm2)
    end
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
  xav = sum(x)/n
  yav = sum(y)/n
  @. x = x - xav
  @. y = y - yav

  # Computing correlation

  tcorr = zeros(Float64,n)
  xnorm = 0.
  ynorm = 0.
  @inbounds for dt in n-1:-1:0
    xnorm = xnorm + x[n-dt]^2
    ynorm = ynorm + y[dt+1]^2
    if xnorm1 != 0 || xnorm2 != 0
      tcorr[dt+1] = LinearAlgebra.dot(@view(x[1:n-dt]),@view(y[dt+1:n]))
      tcorr[dt+1] = tcorr[dt+1] / sqrt(xnorm*ynorm)
    else
      tcorr[dt+1] = 0.
    end
  end

  lags = collect(0:n-1)
  return lags, tcorr

end
