#
# Computes the density function of a list of values
#

function density(v;nbins=nothing,step=nothing,steptype="absolute",vmin=nothing,vmax=nothing)

  ndata = length(v)

  if vmin == nothing
    vmin = minimum(v)
  end
  if vmax == nothing
    vmax = maximum(v)
  end
  if nbins == nothing
    nbins = 100
  end
  if step == nothing
    step = (vmax - vmin)/nbins
  else
    # By default, the step size is absolute
    if steptype == "relative"
      step = step*(vmax-vmin)/nbins
    elseif steptype != "absolute"
      error(" steptype must be \"relative\" or \"absolute\"")
    end
  end

  x = Vector{Float64}(undef,nbins)
  df = Vector{Float64}(undef,nbins)

  binstep = (vmax - vmin)/nbins
  for i in 1:nbins
    x[i] = vmin + (i-1)*binstep + binstep/2
    nv = 0
    for j in 1:ndata
      if ( v[j] > x[i] - step/2 ) && ( v[j] <= x[i] + step/2 )
        nv = nv + 1
      end
    end
    binsize = min(vmax,x[i]+step/2) - max(vmin,x[i]-step/2)
    df[i] = nv/(binsize*ndata)
  end

  return x, df

end

