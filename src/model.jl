using ITensors
using Statistics: mean
include("lattices/cuprate_lattice.jl")
#include("sites/site_hubbholst.jl")
include("sites/threeband/copper.jl")
include("sites/threeband/oxypx.jl")
include("sites/threeband/oxypy.jl")
include("lattices/visualize_lattice.jl")
include("dmrg_lbo.jl")

struct Parameters
    # Model 
    N::Int # number of unit cells
    Nx::Int # number of unit cells along x
    Ny::Int # number of unti cells along y 
    Nsites::Int # number of sites in the MPS representation 
    yperiodic::Bool
    doping::Real
    init_phonons::Int

    # Hamiltonian
    μ::Real
    εd::Real
    εp::Real
    tpd::Real
    tpp::Real
    Vpd::Real
    Upp::Real
    Udd::Real
    ωB1::Real
    ωA1::Real
    gB1::Real
    gA1::Real

    # phonon modes 
    dim_copper_mode_1::Int
    dim_copper_mode_2::Int
    dim_copper_mode_3::Int
    dim_oxygen_x_mode_1::Int
    dim_oxygen_x_mode_2::Int
    dim_oxygen_x_mode_3::Int
    dim_oxygen_y_mode_1::Int
    dim_oxygen_y_mode_2::Int
    dim_oxygen_y_mode_3::Int

    # DMRG 
    DMRG_numsweeps
    DMRG_noise
    DMRG_maxdim
    DMRG_cutoff

    # LBO  
    DMRG_LBO::Bool
    max_LBO_dim::Int
    min_LBO_dim::Int

    # TEBD 
    mid::Int
    T::Int
    τ::Real
    TEBD_cutoff
    TEBD_maxdim
    TEBD_LBO::Bool
end

struct ThreeBandModel
    sites # Lattice
    mpo::MPO # Hamiltonian 
    gates # TEBD Trotterized gates 
end

struct DMRGResults
    nsweep::Int
    maxdim::Vector{Int}
    cutoff::Vector{Float64}
    noise::Vector{Float64}
    ground_state::MPS
    ground_state_energy::Real
    ground_state_entropy::Real
    optimized_basis
    charge_density
    phonon_density
    spin_density
end

mutable struct EquilibriumCorrelations
    start
    stop
    spin
    charge
    particle
    sSC
    dSC_dxdx
    dSC_dpx
    dSC_dydy
    dSC_pyd
    dSC_pypx
    dSC_py1px2
    pSC_dxdx
end

struct TEBDResults
    entropy
    self_overlap
    corrs
    phi_t
    psi_t
end

struct DMRGResultsMinimal
    ground_state_energy
    ground_state_entropy
    charge_density
    phonon_density
    spin_density
end

struct TEBDResultsMinimal
    entropy
    self_overlap
    corrs
end

## SETTING UP THE MODEL ## 
function parameters(; Nx::Int, Ny::Int,
    yperiodic=false,
    doping=0,
    init_phonons::Int=0, μ=0,
    εd=0,
    εp=3,
    tpd=1,
    tpp=0.5,
    Vpd=0,
    Upp=3,
    Udd=8,
    ωB1::Real, # New parameters for the model 
    ωA1::Real,
    gB1::Real,
    gA1::Real,

    # phonon modes 
    dim_copper_mode_1::Int,
    dim_copper_mode_2::Int,
    dim_copper_mode_3::Int,
    dim_oxygen_x_mode_1::Int,
    dim_oxygen_x_mode_2::Int,
    dim_oxygen_x_mode_3::Int,
    dim_oxygen_y_mode_1::Int,
    dim_oxygen_y_mode_2::Int,
    dim_oxygen_y_mode_3::Int, DMRG_numsweeps::Int=20, DMRG_noise=nothing,
    DMRG_maxdim=nothing, DMRG_cutoff=nothing, DMRG_LBO=false,
    max_LBO_dim=nothing, min_LBO_dim=4,
    T::Int=25, τ::Real=0.1, TEBD_cutoff=1E-14, TEBD_maxdim=400, TEBD_LBO=false)

    if isnothing(ωB1)
        ωB1 = 0.5 * tpd
    end
    if isnothing(ωA1)
        ωA1 = 0.05 * tpd
    end
    if isnothing(gB1)
        gB1 = 0.01 * tpd
    end
    if isnothing(gA1)
        gA1 = 0.05 * tpd
    end

    if TEBD_LBO
        @assert DMRG_LBO
    end

    if isnothing(DMRG_noise)
        DMRG_noise = [1E-6, 1E-7, 1E-8, 1E-9, 0]
    end
    if isnothing(DMRG_maxdim)
        DMRG_maxdim = [20, 40, 100, 200, 400]
    end
    if isnothing(DMRG_cutoff)
        DMRG_cutoff = 1E-10
    end
    if isnothing(max_LBO_dim)
        max_LBO_dim = 12
    end
    N = Nx * Ny # number of unit cells 
    Nsites = 3 * N + 2 * Ny # 3 sites per unit cell + another rung to have equal # pos and neg phase bonds
    mid = ceil(Int, Nsites / 2) # midpoint of the DMRG chain 

    return Parameters(N, Nx, Ny, Nsites, yperiodic, doping, init_phonons,
        μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, ωB1, ωA1, gB1, gA1,
        dim_copper_mode_1, dim_copper_mode_2, dim_copper_mode_3,
        dim_oxygen_x_mode_1, dim_oxygen_x_mode_2, dim_oxygen_x_mode_3,
        dim_oxygen_y_mode_1, dim_oxygen_y_mode_2, dim_oxygen_y_mode_3,
        DMRG_numsweeps, DMRG_noise, DMRG_maxdim, DMRG_cutoff, DMRG_LBO,
        max_LBO_dim, min_LBO_dim, mid, T, τ, TEBD_cutoff, TEBD_maxdim, TEBD_LBO)
end

function are_equal(p1::Parameters, p2::Parameters)
    fnames = fieldnames(typeof(p1))
    for fn in fnames
        @assert getfield(p1, fn) == getfield(p2, fn)
    end
end

function calculate_λ(p::Parameters)
    λ = p.gB1^2 / (4 * p.ωB1)
end

## LATTICE ##

function visualize_lattice(p::Parameters)
    Nx, Ny, yperiodic = p.Nx, p.Ny, p.yperiodic
    if yperiodic
        pd_latt = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
        pp_latt = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
        visualize(pd_latt)
        visualize!(pp_latt)
    end
    # even if we are working w periodic boundary conditions, the alternating plus/minus bonds
    # are hidden along the y direction for the pd lattice bc of the periodic bonds, 
    # so we plot it again on top just for visualization's sake 
    pd_latt = OxygenCopper_lattice(Nx, Ny; yperiodic=false, alternate_sign=true)
    pp_latt = OxygenOxygen_lattice(Nx, Ny; yperiodic=false, alternate_sign=true)
    visualize!(pd_latt)
    visualize!(pp_latt)
end

function make_ampo_cuprates(p::Parameters, sites::Vector{Index{Vector{Pair{QN,Int64}}}})
    # Lattice parameters 
    Nsites, Nx, Ny, yperiodic = p.Nsites, p.Nx, p.Ny, p.yperiodic
    # Hubbard parameters 
    μ, εd, εp, tpd, tpp, Vpd, Upp, Udd = p.μ, p.εd, p.εp, p.tpd, p.tpp, p.Vpd, p.Upp, p.Udd
    # Phonon parameters 
    ωB1, ωA1, gB1, gA1 = p.ωB1, p.ωA1, p.gB1, p.gA1
    dim_copper_mode_1, dim_copper_mode_2, dim_copper_mode_3,
    dim_oxygen_x_mode_1, dim_oxygen_x_mode_2, dim_oxygen_x_mode_3,
    dim_oxygen_y_mode_1, dim_oxygen_y_mode_2, dim_oxygen_y_mode_3 = p.dim_copper_mode_1, p.dim_copper_mode_2, p.dim_copper_mode_3,
    p.dim_oxygen_x_mode_1, p.dim_oxygen_x_mode_2, p.dim_oxygen_x_mode_3,
    p.dim_oxygen_y_mode_1, p.dim_oxygen_y_mode_2, p.dim_oxygen_y_mode_3
    dp_lattice = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    pp_lattice = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    site_labels = make_coefficients(Nx + 1, Ny, "Copper", "Oxygen_px", "Oxygen_py")

    # make the hamiltonian 
    ampo = OpSum()

    # on-site terms 
    μ_coefs = make_coefficients(Nx + 1, Ny, εd - μ, εp - μ, εp - μ)
    U_coefs = make_coefficients(Nx + 1, Ny, Udd, Upp, Upp)
    for (n, site_type) in zip(1:Nsites, site_labels)

        # chemical potential term
        ampo .+= μ_coefs[n], "Ntot", n

        # on-site repulsion
        ampo .+= U_coefs[n], "Nupdn", n

        # phonon mode # 1 
        if site_type == "Copper" && dim_copper_mode_1 > 1
            ampo .+= ωB1, "Nb1", n # on-site d mode 1 
        elseif site_type == "Oxygen_px" && dim_oxygen_mode_1 > 1
            ampo .+= ωB1, "Nb1", n # on-site px mode 1 
            ampo .+= gB1, "Ntot(B1d+B1)", n # on-site px coupling mode 1 
        elseif site_type == "Oxygen_py" && dim_oxygen_mode_1 > 1
            ampo .+= ωB1, "Nb1", n # on-site py mode 1 
            ampo .+= -gB1, "Ntot(B1d+B1)", n # on-site py coupling mode 1 - NOTE THE NEGATIVE!
        end

        # phonon mode # 2
        if site_type == "Copper" && dim_copper_mode_2 > 1
            ampo .+= ωA1, "Nb2", n # on-site d mode 2 
        elseif site_type == "Oxygen_px" && dim_oxygen_mode_2 > 1
            ampo .+= ωA1, "Nb2", n # on-site px mode 2
            ampo .+= gA1, "Ntot(B2d+B2)", n # on-site px coupling mode 2
        elseif site_type == "Oxygen_py" && dim_oxygen_mode_2 > 1
            ampo .+= ωA1, "Nb2", n # on-site py mode 2
            ampo .+= gA1, "Ntot(B2d+B2)", n # on-site py coupling mode 2
        end
    end

    # repulsion copper-oxygen
    for b in dp_lattice
        ampo .+= Vpd, "Nup", b.s1, "Nup", b.s2
        ampo .+= Vpd, "Ndn", b.s1, "Ndn", b.s2
        ampo .+= Vpd, "Nup", b.s1, "Ndn", b.s2
        ampo .+= Vpd, "Ndn", b.s1, "Nup", b.s2
    end

    # copper-oxygen hopping
    for b in dp_lattice
        ampo .-= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # oxygen-oxygen hopping
    for b in pp_lattice
        ampo .-= b.sign * tpp, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpp, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpp, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpp, "Cdagdn", b.s2, "Cdn", b.s1
    end

    return MPO(ampo, sites)
end

function number_of_phonon_modes(p::Parameters, ψ::MPS)
    d_idx = to_site_number(p.Ny, 1, 1, "d")
    py_idx = to_site_number(p.Ny, 1, 1, "py")
    px_idx = to_site_number(p.Ny, 1, 1, "px")
    py_site = siteind(ψ, py_idx)
    d_site = siteind(ψ, d_idx)
    px_site = siteind(ψ, px_idx)
    @show dim(d_site)
    @show dim(py_site)
    @show dim(px_site)
end

function make_ampo_cuprates_2mode(p::Parameters, sites::Vector{Index{Vector{Pair{QN,Int64}}}})
    # Lattice parameters 
    Nsites, Nx, Ny, yperiodic = p.Nsites, p.Nx, p.Ny, p.yperiodic
    # Hubbard parameters 
    μ, εd, εp, tpd, tpp, Vpd, Upp, Udd = p.μ, p.εd, p.εp, p.tpd, p.tpp, p.Vpd, p.Upp, p.Udd
    # Phonon parameters 
    ωB1, ωA1, gB1, gA1 = p.ωB1, p.ωA1, p.gB1, p.gA1
    dim_oxygen_x_mode_1, dim_oxygen_y_mode_1 = p.dim_oxygen_x_mode_1, p.dim_oxygen_y_mode_1
    # Lattice
    dp_lattice = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    pp_lattice = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    site_labels = make_coefficients(Nx + 1, Ny, "Copper", "Oxygen_px", "Oxygen_py")

    # make the hamiltonian 
    ampo = OpSum()

    # on-site terms 
    μ_coefs = make_coefficients(Nx + 1, Ny, εd - μ, εp - μ, εp - μ)
    U_coefs = make_coefficients(Nx + 1, Ny, Udd, Upp, Upp)
    for (n, site_type) in zip(1:Nsites, site_labels)

        # chemical potential term
        ampo .+= μ_coefs[n], "Ntot", n

        # on-site repulsion
        ampo .+= U_coefs[n], "Nupdn", n

        # In this set-up, the A1 mode will be on px and the B1 mode will be on py 
        # this is arbitrary!! It doesn't actually matter at all 
        if site_type == "Oxygen_px" && dim_oxygen_x_mode_1 > 1
            ampo .+= ωA1, "Nb1", n # Einstein mode 1 
        elseif site_type == "Oxygen_py" && dim_oxygen_y_mode_1 > 1
            ampo .+= ωB1, "Nb1", n # Einstein mode 2 
        end
    end

    # iterate across unit cells 
    for n in 1:p.N
        # if we have gA1 coupling 
        if gA1 > 0 && dim_oxygen_x_mode_1 > 1
            ph_site = to_site_number(p.Ny, n, "px")
            ex_site = to_site_number(p.Ny, n, "px")
            ey_site = to_site_number(p.Ny, n, "py")
            ampo += gA1, "B1dag+B1", ph_site, "Ntot", ex_site
            ampo += gA1, "B1dag+B1", ph_site, "Ntot", ey_site
        end
        # if we have gB1 coupling 
        if gB1 > 0 && dim_oxygen_y_mode_1 > 1
            ph_site = to_site_number(p.Ny, n, "py")
            ex_site = to_site_number(p.Ny, n, "px")
            ey_site = to_site_number(p.Ny, n, "py")
            ampo += gB1, "B1dag+B1", ph_site, "Ntot", ex_site
            ampo += -gB1, "B1dag+B1", ph_site, "Ntot", ey_site
        end
    end

    # Deal with the end sites also 
    if gB1 > 0 && dim_oxygen_y_mode_1 > 1
        for n in (3*p.N+2):2:p.Nsites
            # n is the (MPS) index of the py site 
            ampo += -gB1, "B1dag+B1", n, "Ntot", n
            @show n
        end
    end

    # repulsion copper-oxygen
    for b in dp_lattice
        ampo .+= Vpd, "Nup", b.s1, "Nup", b.s2
        ampo .+= Vpd, "Ndn", b.s1, "Ndn", b.s2
        ampo .+= Vpd, "Nup", b.s1, "Ndn", b.s2
        ampo .+= Vpd, "Ndn", b.s1, "Nup", b.s2
    end

    # copper-oxygen hopping
    for b in dp_lattice
        ampo .-= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # oxygen-oxygen hopping
    for b in pp_lattice
        ampo .-= b.sign * tpp, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpp, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpp, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpp, "Cdagdn", b.s2, "Cdn", b.s1
    end

    return MPO(ampo, sites)
end

"""
This is the most generic set-up for making a 3BH model  
"""
function make_ampo(p::Parameters, sites::Vector{Index{Vector{Pair{QN,Int64}}}})
    Nsites, Nx, Ny, yperiodic = p.Nsites, p.Nx, p.Ny, p.yperiodic
    μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, ω, g0pp, g0dd, g1pd, g1dp, g1pp = p.μ, p.εd, p.εp, p.tpd, p.tpp, p.Vpd, p.Upp, p.Udd, p.ω, p.g0pp, p.g0dd, p.g1pd, p.g1dp, p.g1pp
    dim_copper_mode_1, dim_copper_mode_2, dim_copper_mode_3,
    dim_oxygen_mode_1, dim_oxygen_mode_2, dim_oxygen_mode_3 = p.dim_copper_mode_1, p.dim_copper_mode_2, p.dim_copper_mode_3,
    p.dim_oxygen_mode_1, p.dim_oxygen_mode_2, p.dim_oxygen_mode_3
    dp_lattice = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    pp_lattice = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    site_labels = make_coefficients(Nx + 1, Ny, "Copper", "Oxygen", "Oxygen")

    # make the hamiltonian 
    ampo = OpSum()

    # on-site terms 
    μ_coefs = make_coefficients(Nx + 1, Ny, εd - μ, εp - μ, εp - μ)
    U_coefs = make_coefficients(Nx + 1, Ny, Udd, Upp, Upp)
    eph_coefs = make_coefficients(Nx + 1, Ny, g0dd, g0pp, g0pp)
    for (n, site_type) in zip(1:Nsites, site_labels)

        # chemical potential term --> INCLUDE 
        ampo .+= μ_coefs[n], "Ntot", n

        # on-site repulsion --> INCLUDE 
        ampo .+= U_coefs[n], "Nupdn", n

        # phonon mode # 1 --> INCLUDE 
        if site_type == "Copper" && dim_copper_mode_1 > 1
            ampo .+= ω, "Nb1", n # Einstein mode 
            ampo .+= eph_coefs[n], "Ntot(B1d+B1)", n # On-site e-ph coupling 
        elseif site_type == "Oxygen" && dim_oxygen_mode_1 > 1
            ampo .+= ω, "Nb1", n
            ampo .+= eph_coefs[n], "Ntot(B1d+B1)", n
        end

        # phonon mode # 2 --> INCLUDE 
        if site_type == "Copper" && dim_copper_mode_2 > 1
            ampo .+= ω, "Nb2", n
            ampo .+= eph_coefs[n], "Ntot(B2d+B2)", n
        elseif site_type == "Oxygen" && dim_oxygen_mode_2 > 1
            ampo .+= ω, "Nb2", n
            ampo .+= eph_coefs[n], "Ntot(B2d+B2)", n
        end

        # phonon mode # 3 --> DONT INCLUDE  
        if site_type == "Copper" && dim_copper_mode_3 > 1
            ampo .+= ω, "Nb3", n
            ampo .+= eph_coefs[n], "Ntot(B3d+B3)", n
        elseif site_type == "Oxygen" && dim_oxygen_mode_3 > 1
            ampo .+= ω, "Nb3", n
            ampo .+= eph_coefs[n], "Ntot(B3d+B3)", n
        end
    end

    # repulsion copper-oxygen --> INCLUDE
    for b in dp_lattice
        ampo .+= Vpd, "Nup", b.s1, "Nup", b.s2
        ampo .+= Vpd, "Ndn", b.s1, "Ndn", b.s2
        ampo .+= Vpd, "Nup", b.s1, "Ndn", b.s2
        ampo .+= Vpd, "Ndn", b.s1, "Nup", b.s2
    end

    # copper-oxygen hopping --> INCLUDE
    for b in dp_lattice
        ampo .-= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # oxygen-oxygen hopping --> INCLUDE
    for b in pp_lattice
        ampo .-= b.sign * tpp, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpp, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpp, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpp, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # electron-phonon nearest neighbour, oxygen-copper 
    for b in dp_lattice
        if site_labels[b.s1] == "Copper"
            # Put ph on copper, put e on oxygen 
            if dim_copper_mode_1 > 1 # --> DONT INCLUDE
                ampo .+= g1dp, "B1dag+B1", b.s1, "Ntot", b.s2
            end
            if dim_copper_mode_2 > 1 # --> DONT INCLUDE
                ampo .+= g1dp, "B2dag+B2", b.s1, "Ntot", b.s2
            end
            if dim_copper_mode_3 > 1 # --> DONT INCLUDE
                ampo .+= g1dp, "B3dag+B3", b.s1, "Ntot", b.s2
            end
            # Put e on copper, put ph on oxygen 
            if dim_oxygen_mode_1 > 1 # --> INCLUDE
                ampo .+= g1pd, "Ntot", b.s1, "B1dag+B1", b.s2
            end
            if dim_oxygen_mode_2 > 1 # --> INCLUDE
                ampo .+= g1pd, "Ntot", b.s1, "B2dag+B2", b.s2
            end
            if dim_oxygen_mode_3 > 1
                ampo .+= g1pd, "Ntot", b.s1, "B3dag+B3", b.s2
            end
        elseif site_labels[b.s1] == "Oxygen"
            # Put ph on oxygen, put e on copper 
            if dim_oxygen_mode_1 > 1
                ampo .+= g1pd, "B1dag+B1", b.s1, "Ntot", b.s2
            end
            if dim_oxygen_mode_2 > 1
                ampo .+= g1pd, "B2dag+B2", b.s1, "Ntot", b.s2
            end
            if dim_oxygen_mode_3 > 1
                ampo .+= g1pd, "B3dag+B3", b.s1, "Ntot", b.s2
            end
            # Put e on oxygen, put ph on copper 
            if dim_copper_mode_1 > 1
                ampo .+= g1dp, "Ntot", b.s1, "B1dag+B1", b.s2
            end
            if dim_copper_mode_2 > 1
                ampo .+= g1dp, "Ntot", b.s1, "B2dag+B2", b.s2
            end
            if dim_copper_mode_3 > 1
                ampo .+= g1dp, "Ntot", b.s1, "B3dag+B3", b.s2
            end
        end
    end

    # copper-oxygen hopping 
    for b in dp_lattice
        ampo .-= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # electron-phonon nearest neighbour, oxygen-oxygen 
    # IS THIS HERMITIAN ??? ## 
    for b in pp_lattice
        # phonon mode # 1
        if dim_oxygen_mode_1 > 1
            # phonon mode # 1 on oxygen, electron on oxygen
            ampo .+= g1pp, "B1dag+B1", b.s1, "Ntot", b.s2
            ampo .+= g1pp, "Ntot", b.s1, "B1dag+B1", b.s2
        end
        if dim_oxygen_mode_2 > 1
            # phonon mode # 2 on oxygen, electron on oxygen
            ampo .+= g1pp, "B2dag+B2", b.s1, "Ntot", b.s2
            ampo .+= g1pp, "Ntot", b.s1, "B2dag+B2", b.s2
        end
        if dim_oxygen_mode_3 > 1
            # phonon mode # 3 on oxygen, electron on oxygen
            ampo .+= g1pp, "B3dag+B3", b.s1, "Ntot", b.s2
            ampo .+= g1pp, "Ntot", b.s1, "B3dag+B3", b.s2
        end
    end

    ## DEFORMATION COUPLING ##
    for b in dp_lattice
        if site_labels[b.s1] == "Copper"
            # Mode 1 (only add phonons to oxygen)
            if dim_oxygen_mode_1 > 1
                ampo .+= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2, "B1dag+B1", b.s2
                ampo .+= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1, "B1dag+B1", b.s2
                ampo .+= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B1dag+B1", b.s2
                ampo .+= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B1dag+B1", b.s2
            end

            # Mode 2 
            if dim_oxygen_mode_2 > 1
                ampo .+= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2, "B2dag+B2", b.s2
                ampo .+= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1, "B2dag+B2", b.s2
                ampo .+= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B2dag+B2", b.s2
                ampo .+= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B2dag+B2", b.s2
            end

            # Mode 3 
            if dim_oxygen_mode_3 > 1
                ampo .+= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2, "B3dag+B3", b.s2
                ampo .+= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1, "B3dag+B3", b.s2
                ampo .+= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B3dag+B3", b.s2
                ampo .+= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B3dag+B3", b.s2
            end

        elseif site_labels[b.s1] == "Oxygen"
            # Mode 1 (only add phonons to oxygen)
            if dim_oxygen_mode_1 > 1
                ampo .+= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2, "B1dag+B1", b.s1
                ampo .+= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1, "B1dag+B1", b.s1
                ampo .+= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B1dag+B1", b.s1
                ampo .+= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B1dag+B1", b.s1
            end

            # Mode 2 
            if dim_oxygen_mode_2 > 1
                ampo .+= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2, "B2dag+B2", b.s1
                ampo .+= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1, "B2dag+B2", b.s1
                ampo .+= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B2dag+B2", b.s1
                ampo .+= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B2dag+B2", b.s1
            end

            # Mode 3 
            if dim_oxygen_mode_3 > 1
                ampo .+= b.sign * tpd, "Cdagup", b.s1, "Cup", b.s2, "B3dag+B3", b.s1
                ampo .+= b.sign * tpd, "Cdagup", b.s2, "Cup", b.s1, "B3dag+B3", b.s1
                ampo .+= b.sign * tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B3dag+B3", b.s1
                ampo .+= b.sign * tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B3dag+B3", b.s1
            end
        end
    end

    return MPO(ampo, sites)
end

"""
Use this function for reconstructing from scratch given site indices 
"""
function ThreeBandModel(p::Parameters, sites::Vector{Index{Vector{Pair{QN,Int64}}}})
    H = make_ampo_cuprates_2mode(p, sites)
    gates = [] #TODO ! make_gates(p, sites)
    # Return the struct 
    ThreeBandModel(sites, H, gates)
end

function ThreeBandModel(p::Parameters, d::DMRGResults)
    sites = siteinds(d.ground_state)
    return ThreeBandModel(p, sites)
end

"""
    ThreeBandModel(p::Parameters)
Make a three-band Hubbard model with phonons given a set of input parameters 
"""
function ThreeBandModel(p::Parameters)
    # sites = siteinds("HubHolst", p.Nsites; dim=p.max_phonons+1)
    # sites = siteinds("TBHSite", p.Nsites)
    site_labels = make_coefficients(Nx + 1, Ny, "Copper", "OxyPx", "OxyPy")
    sites = [siteinds(site_labels[i], 1)[1] for i in 1:p.Nsites]
    return ThreeBandModel(p, sites)
end

## GENERIC STATE FUNCTIONS ##

function compute_overlap(ψ1::MPS, ψ2::MPS)
    LinearAlgebra.norm(inner(ψ1, ψ2))
end

# function compute_phonon_number(ψ::MPS)
#     expect(ψ,"Nb")
# end

function compute_electron_number(ψ::MPS)
    expect(ψ, "Ntot")
end

function compute_entropy(ψ::MPS, b::Int)
    orthogonalize!(ψ, b)
    U, S, V = svd(ψ[b], (linkind(ψ, b - 1), siteind(ψ, b)))
    SvN = 0.0
    for n = 1:dim(S, 1)
        p = S[n, n]^2
        SvN -= p * log(p)
    end
    return SvN
end

## DMRG ## 

function initialize_wavefcn(HM::ThreeBandModel, p::Parameters)
    # Extract the information 
    ups = "Up," * string(p.init_phonons)
    downs = "Dn," * string(p.init_phonons)
    emps = "Emp," * string(p.init_phonons)

    state_arr = make_coefficients(p.Nx + 1, p.Ny, ups, emps, emps)[1:p.Nsites]
    # Make every other hole a down 
    inds_arr = findall(x -> x == ups, state_arr)[begin:2:end]
    state_arr[inds_arr] .= downs
    # Last rung does not get any fermions
    ## SHOULD COMMENT OUT ## 
    state_arr[end-2*Ny+1:end] .= emps

    # Account for doping
    if p.doping > 0
        num_holes_doped = floor(Int, p.doping * p.N)
        @show num_holes_doped
        # put these holes on the px oxygen orbitals, spaced evenly apart 
        spacing = floor(Int, p.N / num_holes_doped) # spacing between unit cells 
        for (i, n) in enumerate(1:spacing:p.N)
            idx = to_site_number(p.Ny, n, "px")
            if isodd(i)
                state_arr[idx] = ups
            else
                state_arr[idx] = downs
            end
        end
    end

    num_holes_up = length(findall(x -> x == ups, state_arr))
    num_holes_down = length(findall(x -> x == downs, state_arr))

    @show num_holes_up
    @show num_holes_down
    #@assert num_holes_up == num_holes_down

    num_holes_total = length(state_arr) - length(findall(x -> x == emps, state_arr))
    @show num_holes_total

    # NOTE: the QN of this state is preserved during DMRG
    productMPS(HM.sites, state_arr)
end

function get_sweeps(p::Parameters, DMRG_numsweeps_per_save::Int)
    # Set sweep params
    sweeps = Sweeps(p.DMRG_numsweeps)
    setnoise!(sweeps, p.DMRG_noise...) # Very important to use noise for this model
    setmaxdim!(sweeps, p.DMRG_maxdim...)
    setmindim!(sweeps, p.DMRG_maxdim...)
    setcutoff!(sweeps, p.DMRG_cutoff...)
    nsweep = p.DMRG_numsweeps
    noise = sweeps.noise
    maxdim = sweeps.maxdim
    cutoff = sweeps.cutoff

    if !isnothing(DMRG_numsweeps_per_save)
        sweeps = Sweeps(DMRG_numsweeps_per_save)
        setnoise!(sweeps, noise[1:DMRG_numsweeps_per_save]...)
        setmaxdim!(sweeps, maxdim[1:DMRG_numsweeps_per_save]...)
        setmindim!(sweeps, maxdim[1:DMRG_numsweeps_per_save]...)
        setcutoff!(sweeps, cutoff[1:DMRG_numsweeps_per_save]...)
        nsweep = DMRG_numsweeps_per_save
        noise = noise[DMRG_numsweeps_per_save:end]
        maxdim = maxdim[DMRG_numsweeps_per_save:end]
        cutoff = cutoff[DMRG_numsweeps_per_save:end]
    end

    return sweeps, nsweep, noise, maxdim, cutoff
end

function get_sweeps(d::DMRGResults, DMRG_numsweeps_per_save::Union{Nothing,Int}=nothing)
    nsweep, noise, maxdim, cutoff = d.nsweep, d.noise, d.maxdim, d.cutoff

    if length(noise) > DMRG_numsweeps_per_save
        if !isnothing(DMRG_numsweeps_per_save)
            sweeps = Sweeps(DMRG_numsweeps_per_save)
            setnoise!(sweeps, noise[1:DMRG_numsweeps_per_save]...)
            setmaxdim!(sweeps, maxdim[1:DMRG_numsweeps_per_save]...)
            setmindim!(sweeps, maxdim[1:DMRG_numsweeps_per_save]...)
            setcutoff!(sweeps, cutoff[1:DMRG_numsweeps_per_save]...)
            nsweep = DMRG_numsweeps_per_save
            noise = noise[DMRG_numsweeps_per_save:end]
            maxdim = maxdim[DMRG_numsweeps_per_save:end]
            cutoff = cutoff[DMRG_numsweeps_per_save:end]
        end
    else
        # Otherwise, we just continue using the same parameters as before, but we augment nsweep
        sweeps = Sweeps(DMRG_numsweeps_per_save)
        setnoise!(sweeps, noise...)
        setmaxdim!(sweeps, maxdim...)
        setmindim!(sweeps, maxdim...)
        setcutoff!(sweeps, cutoff...)
        nsweep += DMRG_numsweeps_per_save
    end

    return sweeps, nsweep, noise, maxdim, cutoff
end

function _run_DMRG(HM::ThreeBandModel, p::Parameters, ϕ0::MPS, sweeps, nsweep, noise, maxdim, cutoff,
    alg, disk_save)

    write_maxdim = maxdim
    if length(maxdim) > 1
        write_maxdim = maxdim[end]
    end

    if p.DMRG_LBO # If performing local basis optimization
        @assert 1 == 0 "Not yet implemented"
        energy, ϕ, Rs = dmrg_lbo(HM.mpo, ϕ0, sweeps, alg=alg, LBO=true,
            max_LBO_dim=p.max_LBO_dim, min_LBO_dim=p.min_LBO_dim)
    else
        if disk_save
            energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg, write_when_maxdim_exceeds=write_maxdim)
        else
            energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg)
        end
        Rs = 0
    end

    # Entropy 
    entropy = compute_entropy(ϕ, ceil(Int, p.Nsites / 2))

    # Electron site densities
    charge_density = ladder_expectation(ϕ, "Ntot", p)
    spin_density = ladder_expectation(ϕ, "Sz", p)

    # Phonon site densities 
    phonon_density = zeros(size(spin_density)..., 9)
    # Copper 
    if p.dim_copper_mode_1 > 1
        phonon_density[:, :, 1] += ladder_expectation(ϕ, "Nb1", p, st="copper")
    end
    if p.dim_copper_mode_2 > 1
        phonon_density[:, :, 2] += ladder_expectation(ϕ, "Nb2", p, st="copper")
    end
    if p.dim_copper_mode_3 > 1
        phonon_density[:, :, 3] += ladder_expectation(ϕ, "Nb3", p, st="copper")
    end
    # Oxygen x
    if p.dim_oxygen_x_mode_1 > 1
        phonon_density[:, :, 4] += ladder_expectation(ϕ, "Nb1", p, st="oxygen_x")
    end
    if p.dim_oxygen_x_mode_2 > 1
        phonon_density[:, :, 5] += ladder_expectation(ϕ, "Nb2", p, st="oxygen_x")
    end
    if p.dim_oxygen_x_mode_3 > 1
        phonon_density[:, :, 6] += ladder_expectation(ϕ, "Nb3", p, st="oxygen_x")
    end
    # Oxygen y
    if p.dim_oxygen_y_mode_1 > 1
        phonon_density[:, :, 7] += ladder_expectation(ϕ, "Nb1", p, st="oxygen_y")
    end
    if p.dim_oxygen_y_mode_2 > 1
        phonon_density[:, :, 8] += ladder_expectation(ϕ, "Nb2", p, st="oxygen_y")
    end
    if p.dim_oxygen_y_mode_3 > 1
        phonon_density[:, :, 9] += ladder_expectation(ϕ, "Nb3", p, st="oxygen_y")
    end

    return DMRGResults(nsweep, maxdim, cutoff, noise, ϕ, energy, entropy,
        Rs, charge_density, phonon_density, spin_density)
end

function run_DMRG(dmrg_results::DMRGResults, HM::ThreeBandModel, p::Parameters;
    DMRG_numsweeps_per_save::Union{Nothing,Int}=nothing,
    alg="divide_and_conquer", disk_save=false, overwrite_sweeps=false)
    # Set DMRG params
    if overwrite_sweeps
        sweeps, nsweep, noise, maxdim, cutoff = get_sweeps(p, DMRG_numsweeps_per_save)
    else
        #start where we left off 
        sweeps, nsweep, noise, maxdim, cutoff = get_sweeps(dmrg_results, DMRG_numsweeps_per_save)
    end

    # Load in the last wavefunction 
    ϕ0 = dmrg_results.ground_state

    _run_DMRG(HM, p, ϕ0, sweeps, nsweep, noise, maxdim, cutoff, alg, disk_save)
end

function run_DMRG(HM::ThreeBandModel, p::Parameters;
    DMRG_numsweeps_per_save::Union{Nothing,Int}=nothing,
    alg="divide_and_conquer", disk_save=false)

    # Set DMRG params
    sweeps, nsweep, noise, maxdim, cutoff = get_sweeps(p, DMRG_numsweeps_per_save)

    # Initialize the wavefunction 
    ϕ0 = initialize_wavefcn(HM, p)

    _run_DMRG(HM, p, ϕ0, sweeps, nsweep, noise, maxdim, cutoff, alg, disk_save)
end

function get_site_indices(st::String, p::Parameters)
    site_labels = make_coefficients(p.Nx + 1, p.Ny, "copper", "oxygen_x", "oxygen_y")[1:p.Nsites]
    indslist = findall(x -> site_labels[x] == st, 1:p.Nsites)
    @assert length(indslist) > 0 "Check that st (sitetype) argument is not capitalized"
    return indslist
end

function ladder_expectation(ϕ::MPS, opname::String, p::Parameters; st::Union{Nothing,String}=nothing)
    density1D = zeros(length(ϕ))
    if isnothing(st)
        indslist = collect(1:length(ϕ))
    else
        indslist = get_site_indices(st, p)
    end
    density1D[indslist] = expect(ϕ, opname, sites=indslist)
    return reshape_into_lattice(density1D, p.Nx, p.Ny)
end

## CORRELATION FUNCTIONS ## 

function apply_onesite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
    ϕ = copy(ϕ)

    ## Account for fermion sign using Jordan-Wigner strings ##
    if opname == "Cup" || opname == "Cdn"
        ϕ = apply_op(ϕ, op(opname, sites[siteidx]), siteidx)
        for i in reverse(1:(siteidx-1)) # Don't act with string on-site
            ϕ = apply_op(ϕ, op("F", sites[i]), i)
        end
        return ϕ

    elseif opname == "Cdagup" || opname == "Cdagdn"
        for i in 1:(siteidx-1) # Don't act with string on-site
            ϕ = apply_op(ϕ, op("F", sites[i]), i)
        end
        ϕ = apply_op(ϕ, op(opname, sites[siteidx]), siteidx)
        return ϕ
    end

    # Otherwise, just apply the operator as usual
    return apply_op(ϕ, op(opname, sites[siteidx]), siteidx)
end

function Π̂(sites, s1::Int, s2::Int)
    # Π = 1/sqrt(2)*(op("Cup",sites[s1])*op("Cdn",sites[s2])
    #                     + op("Cdn",sites[s1])*op("Cup",sites[s2]))

    if s1 > s2 # ordering matters! these are fermions we're talking about
        @warn "s1 > s2"
        s1, s2 = s2, s1
    end

    # First part 
    Π1 = op("Cup", sites[s1])
    for s in (s1+1):(s2-1)
        Π1 *= op("F", sites[s])
    end
    Π1 *= op("Cdn", sites[s2])

    # Second part 
    Π2 = op("Cdn", sites[s1])
    for s in (s1+1):(s2-1)
        Π2 *= op("F", sites[s])
    end
    Π2 *= op("Cup", sites[s2])

    # Take sum and normalize
    Π = 1 / sqrt(2) * (Π1 + Π2)
end

"""
This is my operator 
"""
function Δ̂(sites, s1::Int, s2::Int)
    # Δ_ij = 1/√2 * (c↑i c↓j - c↓i c↑j)

    if s1 > s2 # ordering matters! 
        @warn "s1 > s2"
        s1, s2 = s2, s1
    end

    # First part
    Δ1 = op("Cup", sites[s1])
    for s in (s1+1):(s2-1)
        Δ1 *= op("F", sites[s])
    end
    Δ1 *= op("Cdn", sites[s2])

    # Second part 
    Δ2 = op("Cdn", sites[s1])
    for s in (s1+1):(s2-1)
        Δ2 *= op("F", sites[s])
    end
    Δ2 *= op("Cup", sites[s2])

    # Take difference and normalize
    Δ = 1 / sqrt(2) * (Δ1 - Δ2)
end

function apply_Δ̂(ψ::MPS, sites, s1::Int, s2::Int)
    ψ = copy(ψ)
    if s1 > s2 # ordering matters! 
        @warn "s1 > s2"
        s1, s2 = s2, s1
    end

    # First part
    ψ1 = apply(op("Cup", sites[s1]), ψ)
    for s in (s1+1):(s2-1)
        ψ1 = apply(op("F", sites[s]), ψ1)
    end
    ψ1 = apply(op("Cdn", sites[s2]), ψ1)

    # Second part 
    ψ2 = apply(op("Cdn", sites[s1]), ψ)
    for s in (s1+1):(s2-1)
        ψ2 = apply(op("F", sites[s]), ψ2)
    end
    ψ2 = apply(op("Cup", sites[s2]), ψ2)

    # Take difference and normalize
    Δψ = 1 / sqrt(2) * (ψ1 - ψ2)
end

function apply_Π̂(ψ::MPS, sites, s1::Int, s2::Int)
    # Π = 1/sqrt(2)*(op("Cup",sites[s1])*op("Cdn",sites[s2])
    #                     + op("Cdn",sites[s1])*op("Cup",sites[s2]))

    ψ = copy(ψ)
    if s1 > s2 # ordering matters! these are fermions we're talking about
        @warn "s1 > s2"
        s1, s2 = s2, s1
    end

    # First part 
    ψ1 = apply(op("Cup", sites[s1]), ψ)
    for s in (s1+1):(s2-1)
        ψ1 = apply(op("F", sites[s]), ψ1)
    end
    ψ1 = apply(op("Cdn", sites[s2]), ψ1)

    # Second part 
    ψ2 = apply(op("Cdn", sites[s1]), ψ)
    for s in (s1+1):(s2-1)
        ψ2 = apply(op("F", sites[s]), ψ2)
    end
    ψ2 = apply(op("Cup", sites[s2]), ψ2)

    # Take sum and normalize
    Πψ = 1 / sqrt(2) * (ψ1 + ψ2)
end

"""
This is my code for applying the operator to an MPS at sites s1 and s2 
"""
function apply_twosite_operator(ϕ::MPS, opname::String, sites, s1::Int, s2::Int)
    ϕ = copy(ϕ)
    if opname == "pSC"
        return apply_Π̂(ϕ, sites, s1, s2)
    elseif opname == "dSC"
        return apply_Δ̂(ϕ, sites, s1, s2)
    end
    @error "No recognized two-site operator"
end

function apply_op(ϕ::MPS, op::ITensor, siteidx::Int)
    ϕ = copy(ϕ) # Make a copy of the original state
    orthogonalize!(ϕ, siteidx)
    new_ϕj = op * ϕ[siteidx] # Apply the local operator
    noprime!(new_ϕj)
    ϕ[siteidx] = new_ϕj
    return ϕ
end

function compute_all_equilibrium_correlations(dmrg_results::DMRGResults,
    HM::ThreeBandModel, p::Parameters;
    buffer=nothing)

    if isnothing(buffer)
        buffer = floor(Int, p.Nx / 4)
    end
    # Compute spin and charge correlations
    println("Computing spin correlations for dx-dx bond")
    start, stop, spin_corr = compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "spin", buffer=buffer)

    println("Computing particle-particle correlations for dx-dx bond")
    _, _, charge_corr = compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "particle", buffer=buffer)

    println("Computing charge correlations for dx-dx bond")
    _, _, charge_corr = compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "charge", buffer=buffer)

    # Compute d-wave pairfield correlations 
    pairfield_corrs = compute_equilibrium_pairfield_correlations(dmrg_results, HM, p, buffer=buffer)

    return EquilibriumCorrelations(start, stop, spin_corr, charge_corr, pairfield_corrs...)
end

function compute_equilibrium_pairfield_correlations(dmrg_results::DMRGResults,
    HM::ThreeBandModel, p::Parameters,
    ; buffer=nothing)

    # There are 6 pair-field correlations to compute 
    println("Computing dSC correlations for dx-dx bond")
    _, _, dxdx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dx-dx", "dx-dx", "dSC", buffer=buffer)
    println("Computing dSC correlations for d-px bond")
    _, _, dpx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "d-px", "d-px", "dSC", buffer=buffer)
    println("Computing dSC correlations for dy-dy bond")
    _, _, dydy = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dy-dy", "dy-dy", "dSC", buffer=buffer)
    println("Computing dSC correlations for py-d bond")
    _, _, pyd = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-d", "py-d", "dSC", buffer=buffer)
    println("Computing dSC correlations for py-px bond")
    _, _, pypx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-px", "py-px", "dSC", buffer=buffer)
    println("Computing dSC correlations for py1-px2 bond")
    _, _, py1px2 = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py1-px2", "py1-px2", "dSC", buffer=buffer)

    return [dxdx, dpx, dydy, pyd, pypx, py1px2] # transpose for compatibiity w plotting fcn
end

function get_bonds(bondtype::String, lattice_indices; row=1)
    if bondtype == "d-px"
        s1 = 2
        s2 = 3
        r1, r2 = row, row
    elseif bondtype == "px-d"
        s1 = 3
        s2 = 2
        r1, r2 = row, row
    elseif bondtype == "dx-dx"
        s1 = 2
        s2 = 5
        r1, r2 = row, row
    elseif bondtype == "dy-dy"
        s1 = 2
        s2 = 2
        r1, r2 = 1, 2
    elseif bondtype == "py-d"
        s1 = 1
        s2 = 2
        r1, r2 = row, row
    elseif bondtype == "py1-px2"
        s1 = 1
        s2 = 3
        r1, r2 = 1, 2
    elseif bondtype == "py-px"
        s1 = 1
        s2 = 3
        r1, r2 = row, row
    else
        @error "Bond type not recognized"
        return nothing
    end

    refbond = (lattice_indices[r1, s1], lattice_indices[r2, s2])
    sites_1 = lattice_indices[r1, s1:3:end]
    sites_2 = lattice_indices[r2, s2:3:end]
    bonds = [(sites_1[n], sites_2[n]) for n in 1:length(sites_2)]

    return refbond, bonds
end

function compute_equilibrium_onsite_correlation(dmrg_results::DMRGResults,
    HM::ThreeBandModel, p::Parameters,
    bondtype::String, corrtype::String;
    buffer=nothing, row=1)

    Nx, Ny = p.Nx, p.Ny
    if isnothing(buffer)
        buffer = floor(Int, Nx / 4)
    end
    ϕ = copy(dmrg_results.ground_state)
    Nsites = length(ϕ)
    start = 4 + (buffer - 1) * 3 # discard buffer # of unit cells 
    stop = 3 * Nx - (buffer - 1) * 3 - 2 # discard the last rung and buffer # of unit cells

    # Compute the correlations for each row separately 
    lattice_indices = reshape_into_lattice(collect(1:Nsites), Nx, Ny)
    lattice_indices = convert.(Int, lattice_indices)
    lattice_indices = lattice_indices[:, start:stop]

    # Compute the equilibrium correlations within this chain 
    _, bonds = get_bonds(bondtype, lattice_indices, row=row)
    corr = onsite_correlation(ϕ, bonds, corrtype, HM.sites)

    return start, stop, corr
end

function compute_equilibrium_pairfield_correlation(dmrg_results::DMRGResults,
    HM::ThreeBandModel, p::Parameters,
    bond1::String, bond2::String, SCtype::String;
    buffer=nothing, row=1)

    Nx, Ny = p.Nx, p.Ny
    if isnothing(buffer)
        buffer = floor(Int, Nx / 4)
    end
    ϕ = copy(dmrg_results.ground_state)
    Nsites = length(ϕ)
    start = 4 + (buffer - 1) * 3 # discard buffer # of unit cells 
    stop = 3 * Nx - (buffer - 1) * 3 - 2 # discard the last rung and buffer # of unit cells

    # Compute the correlations for each row separately 
    lattice_indices = reshape_into_lattice(collect(1:Nsites), Nx, Ny)
    lattice_indices = convert.(Int, lattice_indices)
    lattice_indices = lattice_indices[:, start:stop]

    # Compute the equilibrium correlations within this chain 
    refbond, _ = get_bonds(bond1, lattice_indices, row=row)
    _, bonds = get_bonds(bond2, lattice_indices, row=row)
    corr = bond_correlation(ϕ, refbond, bonds, SCtype, HM.sites)

    return start, stop, corr
end

function unzip(bonds)
    sites = []
    for b in bonds
        push!(sites, b[1])
        push!(sites, b[2])
    end
    return unique(sites)
end

function onsite_correlation(ϕ::MPS, bonds, corrtype::String, sites)
    ϕ = copy(ϕ)
    indices = unzip(bonds) # we don't actually care about bonds, bc it's all on single sites
    j = indices[1]

    if corrtype == "spin"
        ψ = apply_onesite_operator(ϕ, "Sz", sites, j)
        corrs = zeros(length(indices))
        Threads.@threads for i in 1:length(indices)
            Szψ = apply_onesite_operator(ψ, "Sz", sites, indices[i])
            corrs[i] = inner(ϕ, Szψ)
        end
        return corrs

    elseif corrtype == "charge"
        ψ = apply_onesite_operator(ϕ, "Ntot", sites, j)
        ninj = zeros(length(indices))
        Threads.@threads for i in 1:length(indices)
            Ntotψ = apply_onesite_operator(ψ, "Ntot", sites, indices[i])
            ninj[i] = inner(ϕ, Ntotψ)
        end
        ni = expect(ϕ, "Ntot")
        nj = ni[j]
        return ninj - nj .* (ni[indices])

    elseif corrtype == "particle"
        ψ = apply_onesite_operator(ϕ, "Cup", sites, j)
        corrs = zeros(length(indices))
        Threads.@threads for i in 1:length(indices)
            cψ = apply_onesite_operator(ψ, "Cup", sites, indices[i])
            corrs[i] = inner(ϕ, cψ)
        end
        return corrs

    elseif corrtype == "sSC"
        ψ = apply_onesite_operator(ϕ, "Cupdn", sites, j)
        corrs = zeros(length(indices))
        Threads.@threads for i in 1:length(indices)
            Σ_iϕ = apply_onesite_operator(ϕ, "Cupdn", sites, indices[i])
            corrs[i] = inner(ψ, Σ_iϕ)
        end
        return corrs
    end

    @error "Unknown correlation type"
end

function bond_correlation(ϕ::MPS, refbond, bonds, corrtype::String, sites)
    if corrtype == "sSC"
        ψ = apply_onesite_operator(ϕ, "Cupdn", sites, j)
        corrs = zeros(length(indices))
        Threads.@threads for i in 1:length(indices)
            Σ_iϕ = apply_onesite_operator(ϕ, "Cupdn", sites, indices[i])
            corrs[i] = inner(ψ, Σ_iϕ)
        end
        return corrs

    elseif corrtype == "pSC"
        ψ = apply_twosite_operator(ϕ, "pSC", sites, refbond[1], refbond[2])
        corrs = zeros(length(bonds))
        Threads.@threads for i in 1:length(bonds)
            Π_iϕ = apply_twosite_operator(ϕ, "pSC", sites, bonds[i][1], bonds[i][2])
            corrs[i] = inner(ψ, Π_iϕ)
        end
        return corrs

    elseif corrtype == "dSC"
        ψ = apply_twosite_operator(ϕ, "dSC", sites, refbond[1], refbond[2])
        corrs = zeros(length(bonds))
        Threads.@threads for i in 1:length(bonds)
            Δ_iϕ = apply_twosite_operator(ϕ, "dSC", sites, bonds[i][1], bonds[i][2])
            corrs[i] = inner(ψ, Δ_iϕ)
        end
        return corrs
    end
    @error "Unknown correlation type"
end