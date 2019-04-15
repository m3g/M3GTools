struct Dihedral

   i :: Int64
   j :: Int64
   k :: Int64
   l :: Int64

   mult :: Int64
   kchi :: Vector{Float64}(undef,4)
   n :: Vector{Int64}(undef,4)
   delta :: Vector{Float64}(undef,4)

end

