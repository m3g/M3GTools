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
 
           end
         end
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
           end
         end
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
           end
         end
       end

       #
       # Read angle parameter data
       #
       # voltar: add dihedrals

     end
     close(file)
   end

end


#  n_bond_lines = 0
#  n_angle_lines = 0
#  n_dihed_lines = 0
#  n_impr_lines = 0
#  n_nonbonded_lines = 0
#  do i = 1, n_parameter_file
#    write(*,*) ' File: ', trim(parameter_file(i))
#    open(10,file=parameter_file(i),action='read',status='old',iostat=ioerr)
#    if ( ioerr /= 0 ) then
#      write(*,*) ' ERROR: Could not open file: ', trim(parameter_file(i))
#      call stop_all()
#    end if
#    read(10,"( a200 )",iostat=ioerr) record
#    if ( ioerr /= 0 ) exit
#    thisfile1 : do
#      if ( comment_line(record) ) then 
#        read(10,"( a200 )",iostat=ioerr) record
#        if ( ioerr /= 0 ) exit thisfile1
#        cycle 
#      end if
#      read(record,*) keyword
#      select case (keyword)
#        case ( "BONDS" ) 
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile1
#              if ( keyword == "ANGLES" .or. &
#                   keyword == "DIHEDRALS" .or. &
#                   keyword == "IMPROPER"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile1
#              end if
#            else
#              n_bond_lines = n_bond_lines + 1
#            end if
#          end do
#        case ( "ANGLES" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, atom3, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile1
#              if ( keyword == "BONDS" .or. &
#                   keyword == "DIHEDRALS" .or. &
#                   keyword == "IMPROPER"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile1
#              end if
#            else
#              n_angle_lines = n_angle_lines + 1
#            end if
#          end do
#        case ( "DIHEDRALS" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, atom3, atom4, dummyr, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile1
#              if ( keyword == "BONDS" .or. &
#                   keyword == "ANGLES" .or. &
#                   keyword == "IMPROPER"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile1
#              end if
#            else
#              n_dihed_lines = n_dihed_lines + 1
#            end if
#          end do
#        case ( "IMPROPER" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, atom3, atom4, dummyr, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile1
#              if ( keyword == "BONDS" .or. &
#                   keyword == "ANGLES" .or. &
#                   keyword == "DIHEDRALS"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile1
#              end if
#            else
#              n_impr_lines = n_impr_lines + 1
#            end if
#          end do
#        case ( "NONBONDED" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, dummyr, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile1
#              if ( keyword == "BONDS" .or. &
#                   keyword == "ANGLES" .or. &
#                   keyword == "DIHEDRALS"  .or. & 
#                   keyword == "IMPROPRER" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile1
#              end if
#            else
#              n_nonbonded_lines = n_nonbonded_lines + 1
#            end if
#          end do
#        case ( "CMAP" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) keyword     
#            if ( ioerr /= 0 ) exit thisfile1
#            if ( keyword == "BONDS" .or. &
#                 keyword == "ANGLES" .or. &
#                 keyword == "DIHEDRALS"  .or. & 
#                 keyword == "IMPROPRER" .or. &
#                 keyword == "NONBONDED" ) then
#              cycle thisfile1
#            end if
#          end do
#        case default
#          read(10,"( a200 )",iostat=ioerr) record
#          if ( ioerr /= 0 ) exit thisfile1
#          cycle thisfile1
#      end select
#    end do thisfile1
#    close(10)
#  end do
#
#  allocate( bond_chars(n_bond_lines,2), &
#            bond_pars(n_bond_lines,2), &
#            angle_chars(n_angle_lines,3), &
#            angle_pars(n_angle_lines,2),&
#            has_urey_bradley(n_angle_lines),&
#            read_urey_bradley(n_angle_lines,2),&
#            dihed_chars(n_dihed_lines,4), &
#            dihed_pars(n_dihed_lines,3), &
#            nonespecific_dihed(ndihed),&
#            impr_chars(n_impr_lines,4), &
#            impr_pars(n_impr_lines,2), &
#            nonespecific_impr(nimpr),&
#            nonbonded_chars(n_nonbonded_lines,1), &
#            nonbonded_pars(n_nonbonded_lines,4) &
#          )
#
#  write(*,*) ' Found ', n_bond_lines, ' bond parameters. '
#  write(*,*) ' Found ', n_angle_lines, ' angle parameters. '
#  write(*,*) ' Found ', n_dihed_lines, ' dihedral parameters. '
#  write(*,*) ' Found ', n_impr_lines, ' improper angle parameters. '
#  write(*,*) ' Found ', n_nonbonded_lines, ' non-bonded parameters. '
#
#  ! Now, the lines containing each type of parameter will be stored in the appropriate
#  ! array
#
#  n_bond_lines = 0
#  n_angle_lines = 0
#  n_dihed_lines = 0
#  n_impr_lines = 0
#  n_nonbonded_lines = 0
#  do i = 1, n_parameter_file
#    open(10,file=parameter_file(i),action='read')
#    read(10,"( a200 )",iostat=ioerr) record
#    thisfile2 : do
#      if ( comment_line(record) ) then 
#        read(10,"( a200 )",iostat=ioerr) record
#        if ( ioerr /= 0 ) exit thisfile2
#        cycle
#      end if
#      read(record,*) keyword
#      select case (keyword)
#        case ( "BONDS" ) 
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile2
#              if ( keyword == "ANGLES" .or. &
#                   keyword == "DIHEDRALS" .or. &
#                   keyword == "IMPROPER"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile2
#              end if
#            else
#              n_bond_lines = n_bond_lines + 1
#              read(record,*) bond_chars(n_bond_lines,1), &
#                             bond_chars(n_bond_lines,2), &
#                             bond_pars(n_bond_lines,1), &
#                             bond_pars(n_bond_lines,2)
#            end if
#          end do
#        case ( "ANGLES" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, atom3, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile2
#              if ( keyword == "BONDS" .or. &
#                   keyword == "DIHEDRALS" .or. &
#                   keyword == "IMPROPER"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile2
#              end if
#            else
#              n_angle_lines = n_angle_lines + 1
#              read(record,*) angle_chars(n_angle_lines,1), &
#                             angle_chars(n_angle_lines,2), &
#                             angle_chars(n_angle_lines,3), &
#                             angle_pars(n_angle_lines,1), &
#                             angle_pars(n_angle_lines,2)
#              ! Check if there are Urey-Bradley terms for this angle
#              has_urey_bradley(n_angle_lines) = .false.
#              read(record,*,iostat=ioerr) (atom1,j=1,3), (dummyr,j=1,4)
#              if ( ioerr == 0 ) then
#                has_urey_bradley(n_angle_lines) = .true.
#                read(record,*,iostat=ioerr) (atom1,j=1,3), (dummyr,j=1,2), &
#                      read_urey_bradley(n_angle_lines,1), &
#                      read_urey_bradley(n_angle_lines,2)
#              end if
#            end if
#          end do
#        case ( "DIHEDRALS" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, atom3, atom4, dummyr, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile2
#              if ( keyword == "BONDS" .or. &
#                   keyword == "ANGLES" .or. &
#                   keyword == "IMPROPER"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile2
#              end if
#            else
#              n_dihed_lines = n_dihed_lines + 1
#              read(record,*) dihed_chars(n_dihed_lines,1), &
#                             dihed_chars(n_dihed_lines,2), &
#                             dihed_chars(n_dihed_lines,3), &
#                             dihed_chars(n_dihed_lines,4), &
#                             dihed_pars(n_dihed_lines,1), &
#                             dihed_pars(n_dihed_lines,2), &
#                             dihed_pars(n_dihed_lines,3)
#            end if
#          end do
#        case ( "IMPROPER" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, atom2, atom3, atom4, dummyr, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile2
#              if ( keyword == "BONDS" .or. &
#                   keyword == "ANGLES" .or. &
#                   keyword == "DIHEDRALS"  .or. & 
#                   keyword == "NONBONDED" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile2
#              end if
#            else
#              n_impr_lines = n_impr_lines + 1
#              read(record,*) impr_chars(n_impr_lines,1), &
#                             impr_chars(n_impr_lines,2), &
#                             impr_chars(n_impr_lines,3), &
#                             impr_chars(n_impr_lines,4), &
#                             impr_pars(n_impr_lines,1), dummyi, &
#                             impr_pars(n_impr_lines,2)
#            end if
#          end do
#        case ( "NONBONDED" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) atom1, dummyr, dummyr, dummyr
#            if ( ioerr /= 0 ) then 
#              read(record,*,iostat=ioerr) keyword     
#              if ( ioerr /= 0 ) exit thisfile2
#              if ( keyword == "BONDS" .or. &
#                   keyword == "ANGLES" .or. &
#                   keyword == "DIHEDRALS"  .or. & 
#                   keyword == "IMPROPRER" .or. &
#                   keyword == "CMAP" ) then
#                cycle thisfile2
#              end if
#            else
#              n_nonbonded_lines = n_nonbonded_lines + 1
#              read(record,*) nonbonded_chars(n_nonbonded_lines,1), dummyr, &
#                             nonbonded_pars(n_nonbonded_lines,1), &
#                             nonbonded_pars(n_nonbonded_lines,2)
#              ! Read 1-4 exclusion parameters if present
#              read(record,*,iostat=ioerr) atom1, (dummyr, j = 1, 6)
#              if ( ioerr == 0 ) then
#                read(record,*) atom1, (dummyr, j = 1, 4), &
#                               nonbonded_pars(n_nonbonded_lines,3), &
#                               nonbonded_pars(n_nonbonded_lines,4)
#              else
#                nonbonded_pars(n_nonbonded_lines,3) = nonbonded_pars(n_nonbonded_lines,1)
#                nonbonded_pars(n_nonbonded_lines,4) = nonbonded_pars(n_nonbonded_lines,2)
#              end if
#            end if
#          end do
#        case ( "CMAP" )
#          do 
#            read(10,"( a200 )",iostat=ioerr) record 
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( comment_line(record) ) cycle
#            read(record,*,iostat=ioerr) keyword     
#            if ( ioerr /= 0 ) exit thisfile2
#            if ( keyword == "BONDS" .or. &
#                 keyword == "ANGLES" .or. &
#                 keyword == "DIHEDRALS"  .or. & 
#                 keyword == "IMPROPRER" .or. &
#                 keyword == "NONBONDED" ) then
#              cycle thisfile2
#            end if
#          end do
#        case default
#          read(10,"( a200 )",iostat=ioerr) record
#          if ( ioerr /= 0 ) exit thisfile2
#          cycle thisfile2
#      end select
#    end do thisfile2
#    close(10)
#  end do
#
#  missing_pars = .false.
#
#  ! Attribute nonbonded parameters to atoms
#
#  write(*,*) ' Setting up non-bonded parameter vectors. '
#  scaling14 = 1.d0
#  atoms: do i = 1, natoms
#    do j = 1, n_nonbonded_lines
#      if ( nonbonded_chars(j,1) == type(i) ) then
#        eps(i) = nonbonded_pars(j,1)
#        sig(i) = nonbonded_pars(j,2)
#        eps14(i) = nonbonded_pars(j,3)
#        sig14(i) = nonbonded_pars(j,4)
#        cycle atoms
#      end if
#    end do
#    write(*,*) ' ERROR: Could not find non-bonded parameters for atom: ', type(i)
#    missing_pars = .true.
#  end do atoms
#
#  ! Attribute bond parameters to bonds
#
#  write(*,*) ' Setting up bond parameter vector. '
#  bonds: do i = 1, nbonds
#    do j = 1, n_bond_lines
#      atom1 = bond_chars(j,1)
#      atom2 = bond_chars(j,2)
#      if ( ( atom1 == type(ibond(i,1)) .and. atom2 == type(ibond(i,2)) ) .or. &
#           ( atom1 == type(ibond(i,2)) .and. atom2 == type(ibond(i,1)) ) ) then
#        fbond(i,1) = bond_pars(j,1)
#        fbond(i,2) = bond_pars(j,2)
#        cycle bonds
#      end if
#    end do
#    write(*,*) ' ERROR: Could not find bond parameters for: ', type(ibond(i,1)), type(ibond(i,2))
#    missing_pars = .true.
#  end do bonds
#
#  ! Attribute angle parameters to angles
#
#  write(*,*) ' Setting up angle parameter vector. '
#
#  ! First, lets count how many Urey-Bradley terms are there
#
#  nureybradley = 0
#  angles1: do i = 1, nangles
#    do j = 1, n_angle_lines
#      atom1 = angle_chars(j,1)
#      atom2 = angle_chars(j,2)
#      atom3 = angle_chars(j,3)
#      if ( ( atom1 == type(iangle(i,1)) .and. &
#             atom2 == type(iangle(i,2)) .and. &
#             atom3 == type(iangle(i,3)) ) .or. &
#           ( atom1 == type(iangle(i,3)) .and. &
#             atom2 == type(iangle(i,2)) .and. & 
#             atom3 == type(iangle(i,1)) ) &
#         ) then
#        if ( has_urey_bradley(j) ) nureybradley = nureybradley + 1
#        cycle angles1
#      end if
#    end do
#  end do angles1
#  allocate( iureybradley(nureybradley,2), fureybradley(nureybradley,2) )
#
#  nureybradley = 0
#  angles2: do i = 1, nangles
#    do j = 1, n_angle_lines
#      atom1 = angle_chars(j,1)
#      atom2 = angle_chars(j,2)
#      atom3 = angle_chars(j,3)
#      if ( ( atom1 == type(iangle(i,1)) .and. &
#             atom2 == type(iangle(i,2)) .and. &
#             atom3 == type(iangle(i,3)) ) .or. &
#           ( atom1 == type(iangle(i,3)) .and. &
#             atom2 == type(iangle(i,2)) .and. & 
#             atom3 == type(iangle(i,1)) ) &
#         ) then
#        fangle(i,1) = angle_pars(j,1) 
#        fangle(i,2) = angle_pars(j,2)*pi/180.d0
#        if ( has_urey_bradley(j) ) then
#          nureybradley = nureybradley + 1
#          iureybradley(nureybradley,1) = iangle(i,1)
#          iureybradley(nureybradley,2) = iangle(i,3)
#          fureybradley(nureybradley,1) = read_urey_bradley(j,1)
#          fureybradley(nureybradley,2) = read_urey_bradley(j,2)
#        end if
#        cycle angles2
#      end if
#    end do
#    write(*,*) ' ERROR: Could not find angle parameters for: ',&
#               type(iangle(i,1)), type(iangle(i,2)), type(iangle(i,3))
#    missing_pars = .true.
#  end do angles2
#
#  ! Dihedral angles: These are read in two steps, the first to allocate the vectors
#  ! which will store the parameters
#
#  ! Summing the number of dihedral parameters per dihedral
#
#  write(*,*) ' Setting up dihedral parameter vector. '
#  do i = 1, ndihed
#    idihed(i,6) = 0
#    nonespecific_dihed(i) = 0
#  end do
#  do j = 1, n_dihed_lines 
#    atom1 = dihed_chars(j,1)
#    atom2 = dihed_chars(j,2)
#    atom3 = dihed_chars(j,3)
#    atom4 = dihed_chars(j,4)
#    do i = 1, ndihed
#      if ( ( atom1 == type(idihed(i,1)) .and. atom2 == type(idihed(i,2)) .and. &
#             atom3 == type(idihed(i,3)) .and. atom4 == type(idihed(i,4)) ) .or. &
#           ( atom1 == type(idihed(i,4)) .and. atom2 == type(idihed(i,3)) .and. &
#             atom3 == type(idihed(i,2)) .and. atom4 == type(idihed(i,1)) ) ) then
#        idihed(i,6) = idihed(i,6) + 1
#      end if
#      if ( ( atom1 == 'X' .and. atom2 == type(idihed(i,2)) .and. &
#             atom3 == type(idihed(i,3)) .and. atom4 == 'X' ) .or. &
#           ( atom1 == 'X' .and. atom2 == type(idihed(i,3)) .and. &
#             atom3 == type(idihed(i,2)) .and. atom4 == 'X' ) ) then 
#        nonespecific_dihed(i) = nonespecific_dihed(i) + 1
#      end if
#    end do
#  end do 
#  do i = 1, ndihed
#    if ( idihed(i,6) == 0 .and. nonespecific_dihed(i) == 0 ) then
#      write(*,*) ' ERROR: Could not find dihedral parameter for: ',&
#                 type(idihed(i,1)), type(idihed(i,2)), type(idihed(i,3)), type(idihed(i,4))
#      missing_pars = .true.
#    end if
#  end do
#
#  ! Setting up the positions in which each dihedral starts in fdihed
#
#  j = 0
#  do i = 1, ndihed
#    idihed(i,5) = j + 1
#    if ( idihed(i,6) == 0 ) then
#      idihed(i,6) = nonespecific_dihed(i)
#    else
#      nonespecific_dihed(i) = 0
#    end if
#    j = j + idihed(i,6)
#  end do
#  nfdihed = j
#
#  ! Allocate fdihed vector ( 3 positions per angle: kchi, n, delta )
#
#  allocate( fdihed(nfdihed,3) )
#
#  ! Fill up fdihed
#
#  k = 0
#  do i = 1, ndihed
#    do j = 1, n_dihed_lines
#      atom1 = dihed_chars(j,1)
#      atom2 = dihed_chars(j,2)
#      atom3 = dihed_chars(j,3)
#      atom4 = dihed_chars(j,4)
#      if ( nonespecific_dihed(i) == 0 ) then 
#        if ( ( atom1 == type(idihed(i,1)) .and. atom2 == type(idihed(i,2)) .and. &
#               atom3 == type(idihed(i,3)) .and. atom4 == type(idihed(i,4)) ) .or. &
#             ( atom1 == type(idihed(i,4)) .and. atom2 == type(idihed(i,3)) .and. &
#               atom3 == type(idihed(i,2)) .and. atom4 == type(idihed(i,1)) ) ) then
#          k = k + 1
#          fdihed(k,1) = dihed_pars(j,1)
#          fdihed(k,2) = dihed_pars(j,2)
#          fdihed(k,3) = dihed_pars(j,3)*pi/180.d0
#        end if
#      else 
#        if ( ( atom1 == 'X' .and. atom2 == type(idihed(i,2)) .and. &
#               atom3 == type(idihed(i,3)) .and. atom4 == 'X' ) .or. &
#             ( atom1 == 'X' .and. atom2 == type(idihed(i,3)) .and. &
#               atom3 == type(idihed(i,2)) .and. atom4 == 'X' ) ) then 
#          k = k + 1
#          fdihed(k,1) = dihed_pars(j,1)
#          fdihed(k,2) = dihed_pars(j,2)
#          fdihed(k,3) = dihed_pars(j,3)*pi/180.d0
#        end if
#      end if
#    end do
#  end do
#
#  ! Improper dihedral angles: These are read in two steps, the first to allocate the vectors
#  ! which will store the parameters, the same as for dihedrals
#
#  ! Summing the number of improper dihedral parameters per dihedral angle
#
#  write(*,*) ' Setting up improper dihedral parameter vector. '
#  do i = 1, nimpr
#    iimpr(i,6) = 0
#    nonespecific_impr(i) = 0
#  end do
#  do j = 1, n_impr_lines 
#    atom1 = impr_chars(j,1)
#    atom2 = impr_chars(j,2)
#    atom3 = impr_chars(j,3)
#    atom4 = impr_chars(j,4)
#    do i = 1, nimpr
#      if ( ( atom1 == type(iimpr(i,1)) .and. atom2 == type(iimpr(i,2)) .and. &
#             atom3 == type(iimpr(i,3)) .and. atom4 == type(iimpr(i,4)) ) .or. &
#           ( atom1 == type(iimpr(i,4)) .and. atom2 == type(iimpr(i,3)) .and. &
#             atom3 == type(iimpr(i,2)) .and. atom4 == type(iimpr(i,1)) ) ) then
#        iimpr(i,6) = iimpr(i,6) + 1
#      end if
#      if ( ( atom1 == type(iimpr(i,1)) .and. atom2 == 'X' .and. &
#             atom3 == 'X' .and. atom4 == type(iimpr(i,4)) ) .or. &
#           ( atom1 == type(iimpr(i,4)) .and. atom2 == 'X' .and. &
#             atom3 == 'X' .and. atom4 == type(iimpr(i,1)) ) ) then 
#        nonespecific_impr(i) = nonespecific_impr(i) + 1
#      end if
#    end do
#  end do 
#
#  ! Setting up the positions in which each improper dihedral starts in fimpr
#
#  j = 0
#  do i = 1, nimpr
#    iimpr(i,5) = j + 1
#    if ( iimpr(i,6) == 0 ) then
#      iimpr(i,6) = nonespecific_impr(i)
#    else
#      nonespecific_impr(i) = 0
#    end if
#    j = j + iimpr(i,6)
#  end do
#  nfimpr = j
#
#  ! Allocate fimpr vector ( 2 positions per angle: kchi, delta )
#
#  allocate( fimpr(nfimpr,2) )
#
#  ! Fill up fimpr
#
#  k = 0
#  do i = 1, nimpr
#    do j = 1, n_impr_lines
#      atom1 = impr_chars(j,1)
#      atom2 = impr_chars(j,2)
#      atom3 = impr_chars(j,3)
#      atom4 = impr_chars(j,4)
#      if ( nonespecific_impr(i) == 0 ) then 
#        if ( ( atom1 == type(iimpr(i,1)) .and. atom2 == type(iimpr(i,2)) .and. &
#               atom3 == type(iimpr(i,3)) .and. atom4 == type(iimpr(i,4)) ) .or. &
#             ( atom1 == type(iimpr(i,4)) .and. atom2 == type(iimpr(i,3)) .and. &
#               atom3 == type(iimpr(i,2)) .and. atom4 == type(iimpr(i,1)) ) ) then
#          k = k + 1
#          fimpr(k,1) = impr_pars(j,1)
#          fimpr(k,2) = impr_pars(j,2) * pi / 180.d0
#        end if
#      else 
#        if ( ( atom1 == type(iimpr(i,1)) .and. atom2 == 'X' .and. &
#               atom3 == 'X' .and. atom4 == type(iimpr(i,4)) ) .or. &
#             ( atom1 == type(iimpr(i,4)) .and. atom2 == 'X' .and.  &
#               atom3 == 'X' .and. atom4 == type(iimpr(i,1)) ) ) then 
#          k = k + 1
#          fimpr(k,1) = impr_pars(j,1)
#          fimpr(k,2) = impr_pars(j,2) * pi / 180.d0
#        end if
#      end if
#    end do
#  end do
#
#  ! Deallocate local vectors
#
#  deallocate( bond_chars, bond_pars, angle_chars, angle_pars, &
#              dihed_chars, dihed_pars, nonespecific_dihed, &
#              impr_chars, impr_pars, nonespecific_impr, &
#              nonbonded_chars, nonbonded_pars, has_urey_bradley, &
#              read_urey_bradley )
#
#  ! Stop if some parameters were missing
# 
#  if ( missing_pars ) call stop_all()
#
#  return
#end subroutine readprm
#
#!
#! Subroutine that tests if a line of a parameter file is an empty
#! or comment line
#!
#
#logical function comment_line(record)
#  integer :: i
#  character(len=200) :: record, record2
#  i = 1
#  do while(record(i:i) <= ' ' .and. i < 200 )
#    i = i + 1
#  end do
#  if ( len_trim(record) == 0 ) then
#    comment_line = .true.
#    return
#  end if
#  record2 = record(i:len(record))
#  comment_line = .false.
#  if ( record2(1:1) == "!" .or. &
#       record2(1:1) == "*" ) then
#    comment_line = .true.
#  end if
#  return
#end function comment_line
          


