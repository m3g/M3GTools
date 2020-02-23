"""

 This function computes the time correlation function of a vector
 defined by the center of masses of two selections.

 The trajectory might be aligned to the first frame according to the
 minimal RMSD of a third selection.

"""

using ProgressMeter

function tcf(simulation :: Simulation,
             abs_start,abs_end,emi_start,emi_end;
             lastframe=0,
             maxdt=-1.,
             align=[],
             theta=-1.,
             r0=1.,
             scaletime=1.)

  # Check if the input is correct for using or not theta

  if theta > 0. 
  
    # This only makes sense if the starting point of the two vectors are the same
  
    if abs_start != emi_start
      error(" Error: If you want to use an angle (theta) to define the emission \n",
            "        vector, the initial point of the emission vector \n",
            "        must be the same as the initial point of the absorption vector. ")
    end

    # And also, this only makes sense if the vectors are not identical

    if ( abs_start == emi_start ) && ( abs_end == emi_end ) 
      error(" Error: If you want to compute the emission vector from a rotation theta \n",
            "        of the absorption vector, the defined absorption and emission \n",
            "        vectors cannot be equal, because they must define the plane \n",
            "        of rotation (actually, this is the only function of the definition \n",
            "        of the emission vector in this case)")
    end
  
    # Converting angle
  
    theta = 3.141592654 * theta / 180. # theta in radians
    cost = cos(theta)
    sint = sin(theta)
    r = Matrix{Float32}(undef,3,3)
  
  end

  nalign = length(align)
  if nalign == 0 
    doalign = false
  else
    doalign = true
    xalign_ref = Matrix{Float32}(undef,nalign+4,3)
    xalign = Matrix{Float32}(undef,nalign+4,3)
  end

  if lastframe == 0 
   lastframe = simulation.nframes
  end
  println(" Last frame to consider = ",lastframe)
  xabs = Matrix{Float32}(undef,lastframe,3)
  xemi = Matrix{Float32}(undef,lastframe,3)

  if maxdt < 0.
    maxstep = lastframe
  else
    maxstep = Int64(maxdt/scaletime)
  end

  # Go to first frame, to be sure
  M3GTools.firstframe(simulation)

  p = Progress(lastframe,5," Reading DCD file: ")
  for iframe in 1:lastframe

    next!(p)
    
    # Reading dcd data for this frame

    sides, xdcd, ydcd, zdcd = M3GTools.nextframe(simulation)

    # Computing the absorption and emission vectors

    cm_abs_start = M3GTools.cm(abs_start,simulation,xdcd,ydcd,zdcd)
    cm_abs_end = M3GTools.cm(abs_end,simulation,xdcd,ydcd,zdcd)
    cm_emi_start = M3GTools.cm(emi_start,simulation,xdcd,ydcd,zdcd)
    cm_emi_end = M3GTools.cm(emi_end,simulation,xdcd,ydcd,zdcd)

    # If align is set, move all vectors according to the alignment of the
    # selected atoms

    if doalign == true 

      if iframe == 1 
        for i in 1:nalign
          xalign_ref[i,1] = xdcd[align[i]]
          xalign_ref[i,2] = ydcd[align[i]]
          xalign_ref[i,3] = zdcd[align[i]]
        end
        for j in 1:3
          xalign_ref[nalign+1,j] = cm_abs_start[j]
          xalign_ref[nalign+2,j] = cm_abs_end[j]
          xalign_ref[nalign+3,j] = cm_emi_start[j]
          xalign_ref[nalign+4,j] = cm_emi_end[j]
        end
      end

      for i in 1:nalign
        xalign[i,1] = xdcd[align[i]]
        xalign[i,2] = ydcd[align[i]]
        xalign[i,3] = zdcd[align[i]]
      end
      for j in 1:3
        xalign[nalign+1,j] = cm_abs_start[j]
        xalign[nalign+2,j] = cm_abs_end[j]
        xalign[nalign+3,j] = cm_emi_start[j]
        xalign[nalign+4,j] = cm_emi_end[j]
      end
 
      xaligned = M3GTools.procrustes(xalign,xalign_ref)

      for j in 1:3
        cm_abs_start[j] = xaligned[nalign+1,j] 
        cm_abs_end[j]   = xaligned[nalign+2,j] 
        cm_emi_start[j] = xaligned[nalign+3,j] 
        cm_emi_end[j]   = xaligned[nalign+4,j] 
      end

    end 

    for i in 1:3
      xabs[iframe,i] = cm_abs_end[i] - cm_abs_start[i]
      xemi[iframe,i] = cm_emi_end[i] - cm_emi_start[i]
    end
    vabs_norm = sqrt( xabs[iframe,1]^2 + xabs[iframe,2]^2 + xabs[iframe,3]^2 )
    vemi_norm = sqrt( xemi[iframe,1]^2 + xemi[iframe,2]^2 + xemi[iframe,3]^2 )
    for i in 1:3
      xabs[iframe,i] = xabs[iframe,i] / vabs_norm
      xemi[iframe,i] = xemi[iframe,i] / vemi_norm
    end

    # Computing emission vector from a rotation of the absoprtion vector, if theta > 0.

    if theta > 0.

      # Rotation axis is perpendicular to both vectors (external product)

      ax = xabs[iframe,2]*xemi[iframe,3] - xabs[iframe,3]*xemi[iframe,2]
      ay = xabs[iframe,3]*xemi[iframe,1] - xabs[iframe,1]*xemi[iframe,3]
      az = xabs[iframe,1]*xemi[iframe,2] - xabs[iframe,2]*xemi[iframe,1]
      vnorm = sqrt( ax^2 + ay^2 + az^2 )
      ax = ax / vnorm
      ay = ay / vnorm
      az = az / vnorm

      # Rotation matrix (quaternion derived)

      xt = 1. - cost
      r[1,1] = cost + ax*ax*xt ; r[1,2] = ax*ay*xt-az*sint; r[1,3] = ax*az*xt+ay*sint
      r[2,1] = ay*ax*xt+az*sint; r[2,2] = cost+ay*ay*xt   ; r[2,3] = ay*az*xt-ax*sint
      r[3,1] = az*ax*xt-ay*sint; r[3,2] = az*ay*xt+ax*sint; r[3,3] = cost + az*az*xt

      # Applying rotation matrix

      xt = xabs[iframe,1]
      yt = xabs[iframe,2]
      zt = xabs[iframe,3]
      xemi[iframe,1] = r[1,1]*xt + r[1,2]*yt + r[1,3]*zt
      xemi[iframe,2] = r[2,1]*xt + r[2,2]*yt + r[2,3]*zt
      xemi[iframe,3] = r[3,1]*xt + r[3,2]*yt + r[3,3]*zt

    end

  end

  return tcf(xabs,xemi,max_dt=maxstep,last_t=lastframe,scaletime=scaletime, r0 = r0, legendre = true)

end

#
# This function computes the time correlation function of two arrays
#

function tcf(v1 :: Array, v2 :: Array; max_dt = -1 ,last_t = -1, scaletime = -1., r0 = 1., legendre = false )

  nv1 = size(v1)[1]
  nv2 = size(v2)[1]
  if length(size(v1)) == 1
    dim = 1
  elseif length(size(v1)) == 2
    dim = size(v1)[2]
    if size(v2)[2] != dim
      error(" For tcf the second and further dimensions of the two vector must be same. ")
    end
  else
    error(" tcf function only accepts two-dimensional arrays. ")
  end

  if last_t == -1
    last_t = min(nv1,nv2)
  end
  if max_dt == -1
    max_dt = last_t - 1
  end
  if scaletime < 0.
    scaletime = 1.
  end

  # Computing the time-dependent correlation function
 
  tcf = zeros(max_dt)
  v1norm = zeros(max_dt)
  v2norm = zeros(max_dt)
  legendre_pol = zeros(max_dt)
  t = Vector{Float32}(undef,max_dt)

  p = Progress(last_t,5," Computing the tcf: ")
  for i in 1:last_t

    next!(p)
    
    v1n = 0.
    for idim in 1:dim
      v1n = v1n + v1[i,idim]*v1[i,idim]
    end

    for j in i:min(i+max_dt-1,last_t)
  
      # Computing the internal product of absoprtion and emission vectors

      int_prod = 0.
      v2n = 0.
      for idim in 1:dim
        int_prod = int_prod + v1[i,idim]*v2[j,idim]
        v2n = v2n + v2[j,idim]*v2[j,idim]
      end

      # Computing the time correlation function and its second legendre polynomial
  
      v1norm[j-i+1] = v1norm[j-i+1] + v1n
      v2norm[j-i+1] = v2norm[j-i+1] + v2n
      tcf[j-i+1] = tcf[j-i+1] + int_prod
      legendre_pol[j-i+1] = legendre_pol[j-i+1] + 0.5*( 3. * int_prod^2 - 1. )
  
    end
  end

  p = Progress(max_dt,5," Final scaling: ")
  for i in 1:max_dt
    next!(p)
    tcf[i] = tcf[i] / ( sqrt(v1norm[i]*v2norm[i]) )# ( last_t - i + 1 )
    legendre_pol[i] = legendre_pol[i] / ( last_t - i + 1 )
    t[i] = scaletime*(i-1)
  end

  if legendre 
    return t, r0*legendre_pol, tcf
  else
    return t, tcf
  end

end

