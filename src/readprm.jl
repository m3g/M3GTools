#
# Subroutine readprm: Reads the parameter files to obtain force constants
#
# L. Martinez, Nov 27, 2013
# Institute of Chemistry - State University of Campinas - Brazil
#
# Tranlsated to Julia in April 15, 2019
#

# Methods to read only one of the components or to call without the flags

readprm!(parfiles, atoms :: Vector{Atom}) = readprm!(parfiles; atoms = atoms) 
readprm!(parfiles, atoms :: Vector{Atom}, bonds :: Vector{Bond}) = readprm!(parfiles; atoms = atoms, bonds = bonds) 
readprm!(parfiles, atoms :: Vector{Atom}, angles :: Vector{Angle}) = readprm!(parfiles; atoms = atoms, angles = angles) 
readprm!(parfiles, atoms :: Vector{Atom}, dihedrals :: Vector{Dihedral}) = readprm!(parfiles; atoms = atoms, dihedrals = dihedrals) 
readprm!(parfiles, atoms :: Vector{Atom}, impropers :: Vector{Improper}) = readprm!(parfiles; atoms = atoms, impropers = impropers) 
readprm!(parfiles, atoms :: Vector{Atom},
                   bonds :: Vector{Bond},
                   angles :: Vector{Angle}, 
                   dihedrals :: Vector{Dihedral},
                   impropers :: Vector{Improper}) = 
  readprm!(parfiles; atoms = atoms, bonds = bonds, angles = angles, dihedrals = dihedrals, impropers = impropers) 

# General function

function readprm!( parfiles; 
                   atoms = nothing, 
                   bonds = nothing, 
                   angles = nothing, 
                   dihedrals = nothing, 
                   impropers = nothing )

   if typeof(parfiles) == String
     parfiles = [ parfiles ]
   end

   dataname = [ "NONBONDED", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "CMAP", "HBOND" ]
   dataread = [    false,    false,    false,     false,       false,     false,  false  ]

   # Vectors for checking if parameters were found
   if atoms != nothing ; check_atoms = zeros(Bool,length(atoms)) ; end
   if bonds != nothing ; check_bonds = zeros(Bool, length(bonds)) ; end
   if angles != nothing ; check_angles = zeros(Bool, length(angles)) ; end
   if dihedrals != nothing ; check_dihedrals = zeros(Bool, length(dihedrals)) ; end

   for file in parfiles

     file = open(file,"r")
     lines = readlines(file)

     for line in lines 

       # Remove comments
       comment = findfirst( isequal('!'), line )
       if comment != nothing
         line = strip(line[1:comment-1])
       else
         line = strip(line)
       end
       if length(line) == 0 ; continue ; end

       # Read data in line
       data = split(line)
       idata = findfirst( isequal(data[1]), dataname )
       if idata != nothing
         @. dataread = false
         dataread[idata] = true
         continue
       end

       #
       # Read atom nonbonded parameter data
       #
       if atoms != nothing && dataread[1]
         for iatom in 1:length(atoms) 
           if data[1] == atoms[iatom].type
             eps = parse(Float64,data[3])
             rmin = parse(Float64,data[4])
             eps14, rmin14 = try
               arse(Float64,data[6]), parse(Float64,data[7])
             catch
                0., 0.
             end
             a = atoms[iatom]
             atoms[iatom] = Atom(a.index,a.residue,a.resid,a.name,a.resname,a.segname,a.type,
                                 a.charge,a.mass,eps,rmin,eps14,rmin14,a.backbone)
             check_atoms[iatom] = true
           end
         end
         continue
       end

       #
       # Read bond parameter data
       #
       if bonds != nothing && dataread[2]
         for ibond in 1:length(bonds) 
           i = bonds[ibond].i
           j = bonds[ibond].j
           if ( ( data[1] == atoms[i].type && data[2] == atoms[j].type ) ||
                ( data[1] == atoms[j].type && data[2] == atoms[i].type ) ) 
             kb = parse(Float64,data[3])
             b0 = parse(Float64,data[4])
             b = bonds[ibond]
             bonds[ibond] = Bond(i,j,kb,b0)
             check_bonds[ibond] = true
           end
         end
         continue
       end

       #
       # Read angle parameter data
       #
       if angles != nothing && dataread[3]
         for iangle in 1:length(angles) 
           i = angles[iangle].i
           j = angles[iangle].j
           k = angles[iangle].k
           if ( ( data[1] == atoms[i].type && data[2] == atoms[j].type && data[3] == atoms[k].type ) ||
                ( data[1] == atoms[k].type && data[2] == atoms[j].type && data[3] == atoms[i].type ) )
             ktheta = parse(Float64,data[4])
             theta0 = parse(Float64,data[5])
             kub, s0 = try
               parse(Float64,data[6]), parse(Float64,data[7])
             catch
               0., 0.
             end
             angles[iangle] = Angle(i,j,k,ktheta,theta0,kub,s0)
             check_angles[iangle] = true
           end
         end
         continue
       end

       #
       # Read dihedral parameter data
       #
       if dihedrals != nothing && dataread[4]
         for idihed in 1:length(dihedrals) 
           i = dihedrals[idihed].i
           j = dihedrals[idihed].j
           k = dihedrals[idihed].k
           l = dihedrals[idihed].l
           if ( data[1] == atoms[i].type && data[2] == atoms[j].type && 
                data[3] == atoms[k].type && data[4] == atoms[l].type ) ||
              ( data[4] == atoms[i].type && data[3] == atoms[j].type && 
                data[2] == atoms[k].type && data[1] == atoms[l].type ) 

             mult = dihedrals[idihed].mult + 1
             kchi = Vector{Float64}(undef,mult)
             n = Vector{Int64}(undef,mult)
             delta = Vector{Float64}(undef,mult)

             if mult > 1 
               kchi[1:dihedrals[idihed].mult] = dihedrals[idihed].kchi
               n[1:dihedrals[idihed].mult] = dihedrals[idihed].n
               delta[1:dihedrals[idihed].mult] = dihedrals[idihed].delta
             end

             kchi[mult] = parse(Float64,data[5])
             n[mult] = parse(Int64,data[6])
             delta[mult] = parse(Float64,data[7])

             dihedrals[idihed] = Dihedral(i,j,k,l,mult,kchi,n,delta)
             check_dihedrals[idihed] = true
           end
         end
         continue
       end

       #
       # Read improper dihedral parameter data
       #
       if impropers != nothing && dataread[5]
         for iimpr in 1:length(impropers) 
           i = impropers[iimpr].i
           j = impropers[iimpr].j
           k = impropers[iimpr].k
           l = impropers[iimpr].l
           if ( ( data[1] == atoms[i].type || data[1] == "X" ) && 
                ( data[2] == atoms[j].type ) && 
                ( data[3] == atoms[k].type ) && 
                ( data[4] == atoms[l].type || data[4] == "X" ) ) ||
              ( ( data[4] == atoms[i].type || data[4] == "X" ) &&
                ( data[3] == atoms[j].type ) && 
                ( data[2] == atoms[k].type ) && 
                ( data[1] == atoms[l].type || data[1] == "X" ) ) ||
              ( ( data[1] == atoms[i].type ) &&
                ( data[2] == atoms[j].type || data[2] == "X" ) && 
                ( data[3] == atoms[k].type || data[3] == "X" ) && 
                ( data[4] == atoms[l].type ) ) ||
              ( ( data[4] == atoms[i].type ) &&
                ( data[3] == atoms[j].type || data[3] == "X" ) && 
                ( data[2] == atoms[k].type || data[2] == "X" ) && 
                ( data[1] == atoms[l].type ) ) 

             mult = impropers[iimpr].mult + 1
             kpsi = Vector{Float64}(undef,mult)
             psi0 = Vector{Int64}(undef,mult)

             if mult > 1 
               kpsi[1:impropers[iimpr].mult] = impropers[iimpr].kpsi
               psi0[1:impropers[iimpr].mult] = impropers[iimpr].psi0
             end

             kpsi[mult] = parse(Float64,data[5])
             psi0[mult] = parse(Float64,data[7])

             impropers[iimpr] = Improper(i,j,k,l,mult,kpsi,psi0)
           end
         end
         continue
       end

     end
     close(file)
   end

   #
   # Report missing parameters
   #
   missing = false
   if atoms != nothing 
     for iatom in 1:length(atoms)
       if ! check_atoms[iatom]
         println(" ERROR: Could not fing NONBONDED parameters for: $(atoms[iatom].type)")
         error = true
       end
     end
   end
   if bonds != nothing 
     for ibond in 1:length(bonds)
       if ! check_bonds[ibond]
         i = bonds[ibond].i
         j = bonds[ibond].j
         println(" ERROR: Could not fing BOND parameters for: $(atoms[i].type) $(atoms[j].type)" )
         missing = true
       end
     end
   end
   if angles != nothing 
     for iangle in 1:length(angles)
       if ! check_angles[iangle]
         i = angles[iangle].i
         j = angles[iangle].j
         k = angles[iangle].k
         println(" ERROR: Could not fing ANGLE parameters for: $(atoms[i].type) $(atoms[j].type) $(atoms[k].type)" )
         missing = true
       end
     end
   end
   if dihedrals != nothing
     for idihed in 1:length(dihedrals)
       if ! check_dihedrals[idihed]
         i = dihedrals[idihed].i
         j = dihedrals[idihed].j
         k = dihedrals[idihed].k
         l = dihedrals[idihed].l
         println(" ERROR: Could not fing DIHEDRAL parameters for: $(atoms[i].type) $(atoms[j].type) $(atoms[k].type) $(atoms[l].type)" )
         missing = true
       end
     end
   end
   if missing
     error(" ERROR: Missing parameters. ")
   end

end

