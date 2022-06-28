
function to_site_number(Ny::Int, x::Int, 
                        y::Int, type::String)
    if type=="d"
        return (x-1)*Ny*3 + (2*y-1)
    elseif type=="py"
        return (x-1)*Ny*3 + 2*y
    elseif type=="px"
        return (x-1)*Ny*3 + 2*Ny + y
    end
end

"""
OxygenCopper_lattice(Nx::Int,
                   Ny::Int;
                   kwargs...)::Lattice
Return a Lattice (array of LatticeBond
objects) corresponding to the two-dimensional
sublattice involving only oxygen-copper bonds 
of dimensions (Nx,Ny).
By default the lattice has periodic boundaries,
but can be made open in the y direction
by specifying the keyword argument
`yperiodic=false`.
"""
function OxygenCopper_lattice(Nx::Int, Ny::Int; kwargs...)::Lattice
  yperiodic = get(kwargs, :yperiodic, true)
  yperiodic = yperiodic && (Ny > 1)
  Nbond = 4*(Nx-2)*(Ny-2) + 4*(Nx-1) + 4*(Ny-2) + 3*(Nx-1) + 3*(Ny-1) + 2 + (yperiodic ? Nx : 0)
  latt = Lattice(undef, Nbond)
  b = 0

  for x in 1:Nx
    for y in 1:Ny
        # The central site 
        n1 = to_site_number(Ny, x, y, "d")

        # Add a pd bond to the right 
        n2 = to_site_number(Ny, x, y, "px")
        latt[b += 1] = LatticeBond(n1, n2, x, y, x, y, "+")

        # Add a pd bond to the left
        if x>1
            n2 = to_site_number(Ny, x-1, y, "px")
            latt[b += 1] = LatticeBond(n1, n2, x, y, x-1, y, "-")
        end

        # Add a pd bond above
        n2 = to_site_number(Ny, x, y, "py")
        latt[b += 1] = LatticeBond(n1, n2, x, y, x, y, "-")

        # Add a pd bond below
        if y>1
            n2 = to_site_number(Ny, x, y-1, "py")
            latt[b += 1] = LatticeBond(n1, n2, x, y, x, y-1, "+")
        elseif yperiodic
            n2 = to_site_number(Ny, x, Ny, "py")
            latt[b += 1] = LatticeBond(n1, n2, x, y, x, Ny, "+")
        end
    end
  end
  return latt
end

"""
OxygenOxygen_lattice(Nx::Int,
                   Ny::Int;
                   kwargs...)::Lattice
Return a Lattice (array of LatticeBond
objects) corresponding to the two-dimensional
sublattice involving only oxygen-oxygen bonds 
of dimensions (Nx,Ny).
By default the lattice has periodic boundaries,
but can be made open in the y direction
by specifying the keyword argument
`yperiodic=false`.
"""
function OxygenOxygen_lattice(Nx::Int, Ny::Int; kwargs...)::Lattice
  yperiodic = get(kwargs, :yperiodic, true)
  yperiodic = yperiodic && (Ny > 1)
  Nbond = 4*(Nx-1)*(Ny-1) + (Ny-1)*2 + (Nx-1)*2 + 1 + (yperiodic ? (Nx-1)*2+1 : 0)
  latt = Lattice(undef, Nbond)
  b = 0

  for x in 1:Nx
    for y in 1:Ny
        # The central site 
        n1 = to_site_number(Ny, x, y, "px")

        # Add a px-py bond to the upper left (in phase)
        n2 = to_site_number(Ny, x, y, "py")
        latt[b += 1] = LatticeBond(n1, n2, x, y, x, y, "+")

        # Add a px-py bond to the upper right (out of phase)
        if x<Nx
            n2 = to_site_number(Ny, x+1, y, "py")
            latt[b += 1] = LatticeBond(n1, n2, x, y, x+1, y, "-")
        end

        # Add a px-py to the lower left (out of phase)
        if y>1
            n2 = to_site_number(Ny, x, y-1, "py")
            latt[b += 1] = LatticeBond(n1, n2, x, y, x, y-1, "-")
        elseif yperiodic
            n2 = to_site_number(Ny, x, Ny, "py")
            latt[b += 1] = LatticeBond(n1, n2, x, y, x, Ny, "-")
        end

        # Add a px-py to the lower right (in phase)
        if x<Nx
            if y>1
                n2 = to_site_number(Ny, x+1, y-1, "py")
                latt[b += 1] = LatticeBond(n1, n2, x, y, x+1, y-1, "+")
            elseif yperiodic
                n2 = to_site_number(Ny, x+1, Ny, "py")
                latt[b += 1] = LatticeBond(n1, n2, x, y, x+1, Ny, "+")
            end
        end
        
    end
  end
  return latt
end