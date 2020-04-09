function sample_interpolation(x,y,xmin,xmax,step)

  nx = length(x)
  n = round(Int64,(xmax-xmin)/step+1)
  xnew = Vector{Float64}(undef,n)
  ynew = Vector{Float64}(undef,n)

  for i in 1:n
   xnew[i] = xmin + (i-1)*step 
   ix = findfirst( x -> x > xnew[i], x )
   if ix == 1 
     ynew[i] = y[1]
   elseif ix == nothing
     ynew[i] = y[nx]
   else
     ynew[i] = y[ix-1] + (xnew[i]-x[ix-1])*((y[ix]-y[ix-1])/(x[ix]-x[ix-1]))
   end
  end

  return xnew, ynew

end
