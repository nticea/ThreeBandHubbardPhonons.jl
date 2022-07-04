"""
    Just a small helper function for converting a coordinates
    in terms of unit cells to a site in the MPS representation
"""
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

struct LatticeBondCuprate
    s1::Int
    s2::Int
    x1::Float64
    y1::Float64
    x2::Float64
    y2::Float64
    sign::Int
    type::String
  end
  
  """
    LatticeBondCuprate(s1::Int,s2::Int,
                    x1::Real,y1::Real,
                    x2::Real,y2::Real,
                    sign::Int,
                    type::String="")
  Construct a LatticeBond struct by
  specifying just the numbers of sites
  1 and 2, or additional details including
  the (x,y) coordinates of the two sites and
  an optional phase factor and/or type string.
  """
  function LatticeBondCuprate(s1::Int, s2::Int)
    return LatticeBondCuprate(s1, s2, 0.0, 0.0, 0.0, 0.0, 1, "")
  end
  
  function LatticeBondCuprate(
    s1::Int, s2::Int, x1::Real, y1::Real, x2::Real, y2::Real, sign::Int, bondtype::String=""
  )
    cf(x) = convert(Float64, x)
    return LatticeBondCuprate(s1, s2, cf(x1), cf(y1), cf(x2), cf(y2), sign, bondtype)
  end
  
  """
  LatticeCuprate is an alias for Vector{LatticeBondCuprate}
  """
  const LatticeCuprate = Vector{LatticeBondCuprate}

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
function OxygenCopper_lattice(Nx::Int, Ny::Int; kwargs...)::LatticeCuprate
  yperiodic = get(kwargs, :yperiodic, true)
  yperiodic = yperiodic && (Ny > 1)
  Nbond = 4*(Nx-2)*(Ny-2) + 4*(Nx-1) + 4*(Ny-2) + 3*(Nx-1) + 3*(Ny-1) + 2 + (yperiodic ? Nx : 0)
  latt = LatticeCuprate(undef, Nbond)
  b = 0

  for x in 1:Nx
    for y in 1:Ny
        # The central site 
        n1 = to_site_number(Ny, x, y, "d")

        # Add a pd bond to the right 
        n2 = to_site_number(Ny, x, y, "px")
        latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, y, 1, "pd")

        # Add a pd bond to the left
        if x>1
            n2 = to_site_number(Ny, x-1, y, "px")
            latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x-1, y, -1, "pd")
        end

        # Add a pd bond above
        n2 = to_site_number(Ny, x, y, "py")
        latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, y, -1, "pd")

        # Add a pd bond below
        if y>1
            n2 = to_site_number(Ny, x, y-1, "py")
            latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, y-1, 1, "pd")
        elseif yperiodic
            n2 = to_site_number(Ny, x, Ny, "py")
            latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, Ny, 1, "pd")
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
function OxygenOxygen_lattice(Nx::Int, Ny::Int; kwargs...)::LatticeCuprate
  yperiodic = get(kwargs, :yperiodic, true)
  yperiodic = yperiodic && (Ny > 1)
  Nbond = 4*(Nx-1)*(Ny-1) + (Ny-1)*2 + (Nx-1)*2 + 1 + (yperiodic ? (Nx-1)*2+1 : 0)
  latt = LatticeCuprate(undef, Nbond)
  b = 0

  for x in 1:Nx
    for y in 1:Ny
        # The central site 
        n1 = to_site_number(Ny, x, y, "px")

        # Add a px-py bond to the upper left (in phase)
        n2 = to_site_number(Ny, x, y, "py")
        latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, y, 1, "pp")

        # Add a px-py bond to the upper right (out of phase)
        if x<Nx
            n2 = to_site_number(Ny, x+1, y, "py")
            latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x+1, y, -1, "pp")
        end

        # Add a px-py to the lower left (out of phase)
        if y>1
            n2 = to_site_number(Ny, x, y-1, "py")
            latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, y-1, -1, "pp")
        elseif yperiodic
            n2 = to_site_number(Ny, x, Ny, "py")
            latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x, Ny, -1, "pp")
        end

        # Add a px-py to the lower right (in phase)
        if x<Nx
            if y>1
                n2 = to_site_number(Ny, x+1, y-1, "py")
                latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x+1, y-1, 1, "pp")
            elseif yperiodic
                n2 = to_site_number(Ny, x+1, Ny, "py")
                latt[b += 1] = LatticeBondCuprate(n1, n2, x, y, x+1, Ny, 1, "pp")
            end
        end
        
    end
  end
  return latt
end

"""
    Make the coefficients based on the positioning of d, px, and py orbitals
"""
function make_coefficients(Nx::Int, Ny::Int, d, px, py)
    coefs = []
    for col in 1:Nx
        for row in 1:Ny
            push!(coefs, d)
            push!(coefs,py)
        end
        for row in 1:Ny
            push!(coefs,px)
        end 
    end
    return coefs 
end

"""
    function argsort_lattice(a, Nx::Int, Ny::Int)
given an array a of length Nsites, it will rearrange it into the shape of the lattice
the unit cells are flattened along the x direction as (py , d , px)
returns an array of shape (Ny x 3*Nx)
"""
function reshape_into_lattice(a, Nx::Int, Ny::Int)
    a_lattice = zeros(Ny, 3*Nx)
    # Iterate through each of the cells 
    for y in 1:Ny
        s=0
        for x in 1:Nx
            idx_d = to_site_number(Ny, x, y, "d")
            idx_px = to_site_number(Ny, x, y, "px")
            idx_py = to_site_number(Ny, x, y, "py")
            a_lattice[y, s+=1] = a[idx_py]
            a_lattice[y, s+=1] = a[idx_d]
            a_lattice[y, s+=1] = a[idx_px]
        end
    end
    return a_lattice 
end
