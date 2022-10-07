using LinearAlgebra
using SparseArrays

"""
    space(::SiteType"TBHSite"; 
          conserve_qns = false,
          conserve_sz = conserve_qns,
          conserve_nf = conserve_qns,
          conserve_nfparity = conserve_qns,
          qnname_sz = "Sz",
          qnname_nf = "Nf",
          qnname_nfparity = "NfParity")

Create the Hilbert space for a site of type "HubHolst".

Optionally specify the conserved symmetries and their quantum number labels.
"""

## GLOBAL CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
DIM_1 = 2 # maximum 1 phonon 
DIM_2 = 2 # maximum 1 phonon 
DIM_3 = 2 # maximum 1 phonon 

function ITensors.space(
  ::SiteType"TBHSite";
  conserve_qns=true,
  conserve_sz=conserve_qns,
  conserve_nf=conserve_qns,
  conserve_nfparity=conserve_qns,
  qnname_sz="Sz",
  qnname_nf="Nf",
  qnname_nfparity="NfParity")

  dim = DIM_1 * DIM_2 * DIM_3 

  if conserve_sz && conserve_nf
    return [
      QN((qnname_nf, 0, -1), (qnname_sz, 0)) => 1 * dim 
      QN((qnname_nf, 1, -1), (qnname_sz, +1)) => 1 * dim
      QN((qnname_nf, 1, -1), (qnname_sz, -1)) => 1 * dim
      QN((qnname_nf, 2, -1), (qnname_sz, 0)) => 1 * dim
    ]
  elseif conserve_nf
    return [
      QN(qnname_nf, 0, -1) => 1 * dim
      QN(qnname_nf, 1, -1) => 2 * dim
      QN(qnname_nf, 2, -1) => 1 * dim
    ]
  elseif conserve_sz
    return [
      QN((qnname_sz, 0), (qnname_nfparity, 0, -2)) => 1 * dim
      QN((qnname_sz, +1), (qnname_nfparity, 1, -2)) => 1 * dim
      QN((qnname_sz, -1), (qnname_nfparity, 1, -2)) => 1 * dim
      QN((qnname_sz, 0), (qnname_nfparity, 0, -2)) => 1 * dim
    ]
  elseif conserve_nfparity
    return [
      QN(qnname_nfparity, 0, -2) => 1 * dim
      QN(qnname_nfparity, 1, -2) => 2 * dim
      QN(qnname_nfparity, 0, -2) => 1 * dim
    ]
  end
  return 4 * dim 
end

function ITensors.val(::ValName{N}, ::SiteType"TBHSite") where {N}
    hubbard_type, phonon_num = split(String(N), ",")
    if hubbard_type == "Emp"
        n1 = 1
    elseif hubbard_type == "Up"
        n1 = 2
    elseif hubbard_type == "Dn"
        n1 = 3
    elseif hubbard_type == "UpDn"
        n1 = 4
    else
        throw(DomainError(hubbard_type, "expects Emp, Up, Dn, UpDn"))
    end

    n2 = parse(Int, String(phonon_num)) + 1
    return n1+n2
end

function ITensors.state(n::StateName{N}, ::SiteType"TBHSite", s::Index) where {N}
    hubbard_type, n = split(String(ITensors.name(n)), ",")
    n = parse(Int, n)

    if hubbard_type == "Emp"
        st_hubb = [1, 0, 0, 0]
    elseif hubbard_type == "Up"
        st_hubb = [0, 1, 0, 0]
    elseif hubbard_type == "Dn"
        st_hubb = [0, 0, 1, 0]
    elseif hubbard_type == "UpDn"
        st_hubb = [0, 0, 0, 1]
    else
        throw(DomainError(hubbard_type, "expects Emp, Up, Dn, UpDn"))
    end

    st_ph = zeros(floor(Int,dim(s)/4))
    st_ph[n+1] = 1

    return kron(st_hubb, st_ph)
end

## ELECTRON OPERATORS ## 

function ITensors.op(::OpName"I", ::SiteType"TBHSite", s::Index)
  M = Matrix(I, Int(dim(s)), Int(dim(s)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nup", ::SiteType"TBHSite", s::Index)
  Nup = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
  ]
  M = kron(Nup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ndn", ::SiteType"TBHSite", s::Index)
  Ndn = [
      0.0 0.0 0.0 0.0
      0.0 0.0 0.0 0.0
      0.0 0.0 1.0 0.0
      0.0 0.0 0.0 1.0
    ]
  M = kron(Ndn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nupdn", ::SiteType"TBHSite", s::Index)
  Nupdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
  ]
  M = kron(Nupdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ntot", ::SiteType"TBHSite", s::Index)
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ]
  M = kron(Ntot, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cup", ::SiteType"TBHSite", s::Index)
  Cup = [
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Cup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdagup", ::SiteType"TBHSite", s::Index)
  Cdagup = [
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
  ]
  M = kron(Cdagup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdn", ::SiteType"TBHSite", s::Index)
  Cdn = [
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Cdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cupdn", ::SiteType"TBHSite", s::Index)
  Cupdn = [
    0.0  0.0  0.0  -1.0
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
  ]
  M = kron(Cupdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdagupdn", ::SiteType"TBHSite", s::Index)
  Cdagupdn = [
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
    0.0  0.0  0.0   0.0
   -1.0  0.0  0.0   0.0
  ]
  M = kron(Cdagupdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Cdagdn", ::SiteType"TBHSite", s::Index)
  Cdagdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
  ]
  M = kron(Cdagdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Aup", ::SiteType"TBHSite", s::Index)
  Aup = [
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 1.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Aup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Adagup", ::SiteType"TBHSite", s::Index)
  Adagup = [
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
  ]
  M = kron(Adagup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Adn", ::SiteType"TBHSite", s::Index)
  Adn = [
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Adn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Adagdn", ::SiteType"TBHSite", s::Index)
  Adagdn = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
  ]
  M = kron(Adagdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"F", ::SiteType"TBHSite", s::Index)
  F = [
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
    0.0 0.0 -1.0 0.0
    0.0 0.0 0.0 1.0
  ]
  M = kron(F, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Fup", ::SiteType"TBHSite", s::Index)
  Fup = [
    1.0 0.0 0.0 0.0
    0.0 -1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 -1.0
  ]
  M = kron(Fup, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Fdn", ::SiteType"TBHSite", s::Index)
  Fdn = [
    1.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 -1.0 0.0
    0.0 0.0 0.0 -1.0
  ]
  M = kron(Fdn, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Sz", ::SiteType"TBHSite", s::Index)
  Sz = [
    0.0 0.0 0.0 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 -0.5 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Sz, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Sx", ::SiteType"TBHSite", s::Index)
  Sx = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.5 0.0
    0.0 0.5 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Sx, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"S+", ::SiteType"TBHSite", s::Index)
  Splus = [
    0.0 0.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Splus, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"S-", ::SiteType"TBHSite", s::Index)
  Sminus = [
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 0.0 0.0
  ]
  M = kron(Sminus, Matrix(I, floor(Int,dim(s)/4), floor(Int,dim(s)/4)))
  ITensor(M, prime(s), dag(s))
end

ITensors.has_fermion_string(::OpName"Cup", ::SiteType"TBHSite") = true
function ITensors.has_fermion_string(on::OpName"c↑", st::SiteType"TBHSite")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdagup", ::SiteType"TBHSite") = true
function ITensors.has_fermion_string(on::OpName"c†↑", st::SiteType"TBHSite")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdn", ::SiteType"TBHSite") = true
function ITensors.has_fermion_string(on::OpName"c↓", st::SiteType"TBHSite")
  return has_fermion_string(alias(on), st)
end
ITensors.has_fermion_string(::OpName"Cdagdn", ::SiteType"TBHSite") = true
function ITensors.has_fermion_string(on::OpName"c†↓", st::SiteType"TBHSite")
  return has_fermion_string(alias(on), st)
end

## PHONON OPERATORS ## 

## MODE 1 ## 
function ITensors.op(::OpName"B1", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_1 > 1 "This mode has no phonons"

  # Define the lowering operator
  B = Tridiagonal(zeros(DIM_1-1),zeros(DIM_1),sqrt.(collect(1:DIM_1-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), B, Matrix(I, DIM_2, DIM_2), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"B1dag", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_1 > 1 "This mode has no phonons"

  # Define the raising operator
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_1-1)),zeros(DIM_1),zeros(DIM_1-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Bdag, Matrix(I, DIM_2, DIM_2), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb1", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_1 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_1-1)),zeros(DIM_1),zeros(DIM_1-1)))
  B = Tridiagonal(zeros(DIM_1-1),zeros(DIM_1),sqrt.(collect(1:DIM_1-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), (Bdag * B), Matrix(I, DIM_2, DIM_2), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb1^2", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_1 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_1-1)),zeros(DIM_1),zeros(DIM_1-1)))
  B = Tridiagonal(zeros(DIM_1-1),zeros(DIM_1),sqrt.(collect(1:DIM_1-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), (Bdag * B * Bdag * B), Matrix(I, DIM_2, DIM_2), 
                                                  Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"B1dag+B1", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_1 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_1-1)),zeros(DIM_1),zeros(DIM_1-1)))
  B = Tridiagonal(zeros(DIM_1-1),zeros(DIM_1),sqrt.(collect(1:DIM_1-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), (Bdag + B), Matrix(I, DIM_2, DIM_2), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ntot(B1d+B1)", ::SiteType"TBHSite", s::Index)
  # Fermion part 
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ] 
  
  # Make B†+B
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons" 
  @assert DIM_1 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_1-1)),zeros(DIM_1),zeros(DIM_1-1))) 
  B = Tridiagonal(zeros(DIM_1-1),zeros(DIM_1),sqrt.(collect(1:DIM_1-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Ntot, (Bdag + B), Matrix(I, DIM_2, DIM_2), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

## MODE 2 ## 
function ITensors.op(::OpName"B2", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_2 > 1 "This mode has no phonons"

  # Define the lowering operator
  B = Tridiagonal(zeros(DIM_2-1),zeros(DIM_2),sqrt.(collect(1:DIM_2-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), B, Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"B2dag", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_2 > 1 "This mode has no phonons"

  # Define the raising operator
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_2-1)),zeros(DIM_2),zeros(DIM_2-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), Bdag, Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb2", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_2 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_2-1)),zeros(DIM_2),zeros(DIM_2-1)))
  B = Tridiagonal(zeros(DIM_2-1),zeros(DIM_2),sqrt.(collect(1:DIM_2-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), (Bdag * B), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb2^2", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_2 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_2-1)),zeros(DIM_2),zeros(DIM_2-1)))
  B = Tridiagonal(zeros(DIM_2-1),zeros(DIM_2),sqrt.(collect(1:DIM_2-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), (Bdag * B * Bdag * B), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"B2dag+B2", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_2 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_2-1)),zeros(DIM_2),zeros(DIM_2-1)))
  B = Tridiagonal(zeros(DIM_2-1),zeros(DIM_2),sqrt.(collect(1:DIM_2-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), (Bdag + B), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ntot(B2d+B2)", ::SiteType"TBHSite", s::Index)
  # Fermion part 
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ] 
  
  # Make B†+B
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons" 
  @assert DIM_2 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_2-1)),zeros(DIM_2),zeros(DIM_2-1))) 
  B = Tridiagonal(zeros(DIM_2-1),zeros(DIM_2),sqrt.(collect(1:DIM_2-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Ntot, Matrix(I, DIM_1, DIM_1), (Bdag + B), Matrix(I, DIM_3, DIM_3))
  ITensor(M, prime(s), dag(s))
end

## MODE 3 ## 
function ITensors.op(::OpName"B3", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_3 > 1 "This mode has no phonons"

  # Define the lowering operator
  B = Tridiagonal(zeros(DIM_3-1),zeros(DIM_3),sqrt.(collect(1:DIM_3-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), Matrix(I, DIM_2, DIM_2), B)
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"B3dag", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_3 > 1 "This mode has no phonons"

  # Define the raising operator
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_3-1)),zeros(DIM_3),zeros(DIM_3-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), Matrix(I, DIM_2, DIM_2), Bdag)
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb3", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_3 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_3-1)),zeros(DIM_3),zeros(DIM_3-1)))
  B = Tridiagonal(zeros(DIM_3-1),zeros(DIM_3),sqrt.(collect(1:DIM_3-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), Matrix(I, DIM_2, DIM_2), (Bdag * B))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Nb3^2", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_3 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_3-1)),zeros(DIM_3),zeros(DIM_3-1)))
  B = Tridiagonal(zeros(DIM_3-1),zeros(DIM_3),sqrt.(collect(1:DIM_3-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), Matrix(I, DIM_2, DIM_2), 
                                                      (Bdag * B * Bdag * B))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"B3dag+B3", ::SiteType"TBHSite", s::Index)
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons"
  @assert DIM_3 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_3-1)),zeros(DIM_3),zeros(DIM_3-1)))
  B = Tridiagonal(zeros(DIM_3-1),zeros(DIM_3),sqrt.(collect(1:DIM_3-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Matrix(I, 4, 4), Matrix(I, DIM_1, DIM_1), Matrix(I, DIM_2, DIM_2), (Bdag + B))
  ITensor(M, prime(s), dag(s))
end

function ITensors.op(::OpName"Ntot(B3d+B3)", ::SiteType"TBHSite", s::Index)
  # Fermion part 
  Ntot = [
    0.0 0.0 0.0 0.0
    0.0 1.0 0.0 0.0
    0.0 0.0 1.0 0.0
    0.0 0.0 0.0 2.0
  ] 
  
  # Make B†+B
  @assert floor(Int,dim(s)/4) > 1 "Must allow a nonzero number of phonons" 
  @assert DIM_3 > 1 "This mode has no phonons"

  # Define the operators
  Bdag = sparse(Tridiagonal(sqrt.(collect(1:DIM_3-1)),zeros(DIM_3),zeros(DIM_3-1))) 
  B = Tridiagonal(zeros(DIM_3-1),zeros(DIM_3),sqrt.(collect(1:DIM_3-1)))
  
  # Take Kronecker product with all the appropriate subspaces 
  M = kron(Ntot, Matrix(I, DIM_1, DIM_1), Matrix(I, DIM_2, DIM_2), (Bdag + B))
  ITensor(M, prime(s), dag(s))
end