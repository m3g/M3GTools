function distribution(v;nbins=Nothing,step=Nothing,vmin=Nothing,vmax=Nothing,norm=0)

  ndata = length(v)

  if vmin == Nothing
    vmin = minimum(v)
  end
  if vmax == Nothing
    vmax = maximum(v)
  end
  if step == Nothing && nbins == Nothing
    step = (vmax - vmin)/(ndata/50)
    nbins = round(Int64,(vmax-vmin)/step)
  elseif step != Nothing && nbins == Nothing
    nbins = round(Int64,(vmax-vmin)/step)
  elseif step == Nothing && nbins != Nothing
    step = (vmax-vmin)/nbins
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
    df[i] = nv
  end

  # Normalize

  ntotal = 0
  for i in 1:nbins
    ntotal = ntotal + df[i]
  end
  df = df / ntotal

  return x, df

end

