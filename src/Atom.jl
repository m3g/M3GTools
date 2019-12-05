struct Atom

  index :: Int64
  residue :: Int64
  resid :: Int64
  
  name :: String
  resname :: String
  segname :: String
  type :: String
  element :: String
  chain :: String

  charge :: Float64
  mass :: Float64
  eps :: Float64
  rmin :: Float64
  eps14 :: Float64
  rmin14 :: Float64
  
  backbone :: Bool

  coor :: Vector{Float64}

end
