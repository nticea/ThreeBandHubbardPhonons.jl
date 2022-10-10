using ITensors
using Statistics: mean
include("lattices/cuprate_lattice.jl")
#include("sites/site_hubbholst.jl")
include("sites/threeband/copper.jl")
include("sites/threeband/oxygen.jl")
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
    ω::Real
    g0pp::Real
    g0dd::Real
    g1pd::Real
    g1dp::Real
    g1pp::Real

    # phonon modes 
    dim_copper_mode_1::Int
    dim_copper_mode_2::Int
    dim_copper_mode_3::Int
    dim_oxygen_mode_1::Int
    dim_oxygen_mode_2::Int
    dim_oxygen_mode_3::Int

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
    dSC_dxdx
    dSC_dpx
    dSC_dydy
    dSC_pyd
    dSC_pypx
    dSC_py1px2
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
function parameters(;Nx::Int, Ny::Int, 
    yperiodic=false, 
    doping=0,
    init_phonons::Int=0,

    μ=0,
    εd=0,
    εp=3,
    tpd=1,
    tpp=0.5,
    Vpd=0,
    Upp=3,
    Udd=8,
    ω::Real,
    g0pp::Real,
    g0dd::Real,
    g1pd::Real,
    g1dp::Real,
    g1pp::Real,

    # phonon modes 
    dim_copper_mode_1::Int,
    dim_copper_mode_2::Int,
    dim_copper_mode_3::Int,
    dim_oxygen_mode_1::Int,
    dim_oxygen_mode_2::Int,
    dim_oxygen_mode_3::Int,

    DMRG_numsweeps::Int=20, DMRG_noise=nothing, 
    DMRG_maxdim=nothing, DMRG_cutoff=nothing, DMRG_LBO=false,
    max_LBO_dim=nothing, min_LBO_dim=4,
    T::Int=25, τ::Real=0.1, TEBD_cutoff=1E-14, TEBD_maxdim=400, TEBD_LBO=false)

    if isnothing(ω)
        ω = 0.5*tpd
    end
    if isnothing(g0pp)
        g0pp = 0.05*tpd
    end
    if isnothing(g0dd)
        g0dd = 0.05*tpd
    end
    if isnothing(g1pd)
        g1pd = 0
    end
    if isnothing(g1dp)
        g1dp = 0
    end
    if isnothing(g1pp)
        g1pp = 0
    end

    if TEBD_LBO
        @assert DMRG_LBO
    end

    if isnothing(DMRG_noise)
        DMRG_noise = [1E-6,1E-6,1E-8,0]
    end
    if isnothing(DMRG_maxdim)
        DMRG_maxdim = [20,40,100,200,400]
    end
    if isnothing(DMRG_cutoff)
        DMRG_cutoff = 1E-10
    end
    if isnothing(max_LBO_dim)
        max_LBO_dim = 12
    end
    N = Nx*Ny # number of unit cells 
    Nsites = 3*N + 2*Ny # 3 sites per unit cell + another rung to have equal # pos and neg phase bonds
    mid = ceil(Int,Nsites/2) # midpoint of the DMRG chain 

    return Parameters(N, Nx, Ny, Nsites, yperiodic, doping, init_phonons,
            μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, ω, g0pp, g0dd, g1pd, g1dp, g1pp,
            dim_copper_mode_1, dim_copper_mode_2, dim_copper_mode_3,
            dim_oxygen_mode_1, dim_oxygen_mode_2, dim_oxygen_mode_3,
            DMRG_numsweeps, DMRG_noise, DMRG_maxdim, DMRG_cutoff, DMRG_LBO,
            max_LBO_dim, min_LBO_dim, mid, T, τ, TEBD_cutoff, TEBD_maxdim, TEBD_LBO)
end

function are_equal(p1::Parameters, p2::Parameters)
    fnames = fieldnames(typeof(p1))
    for fn in fnames 
        @assert getfield(p1, fn)==getfield(p2, fn)
    end
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

function make_ampo(p::Parameters, sites::Vector{Index{Vector{Pair{QN, Int64}}}})
    Nsites, Nx, Ny, yperiodic = p.Nsites, p.Nx, p.Ny, p.yperiodic
    μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, ω, g0pp, g0dd, g1pd, g1dp, g1pp = p.μ, p.εd, p.εp, p.tpd, p.tpp, p.Vpd, p.Upp, p.Udd, p.ω, p.g0pp, p.g0dd, p.g1pd, p.g1dp, p.g1pp
    dim_copper_mode_1, dim_copper_mode_2, dim_copper_mode_3,
            dim_oxygen_mode_1, dim_oxygen_mode_2, dim_oxygen_mode_3 = p.dim_copper_mode_1, p.dim_copper_mode_2, p.dim_copper_mode_3,
                                                                    p.dim_oxygen_mode_1, p.dim_oxygen_mode_2, p.dim_oxygen_mode_3
    dp_lattice = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    pp_lattice = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic, alternate_sign=true)
    site_labels = make_coefficients(Nx+1, Ny, "Copper", "Oxygen", "Oxygen")

    # make the hamiltonian 
    ampo = OpSum()

    # on-site terms 
    μ_coefs = make_coefficients(Nx+1, Ny, εd-μ, εp-μ, εp-μ)
    U_coefs = make_coefficients(Nx+1, Ny, Udd, Upp, Upp)
    eph_coefs = make_coefficients(Nx+1, Ny, g0dd, g0pp, g0pp)
    for (n, site_type) in zip(1:Nsites, site_labels)
        
        # chemical potential term
        # NOTE: MAYBE THIS IS DOING UNNECESSARY WORK -- add Δ term only to px/py
        ampo .+= μ_coefs[n], "Ntot", n
        
        # on-site repulsion
        ampo .+= U_coefs[n], "Nupdn", n

        # phonon mode # 1
        if site_type == "Copper" && dim_copper_mode_1 > 1
            ampo .+= ω, "Nb1", n # Einstein mode 
            ampo .+= eph_coefs[n], "Ntot(B1d+B1)", n # On-site e-ph coupling 
        elseif site_type == "Oxygen" && dim_oxygen_mode_1 > 1
            ampo .+= ω, "Nb1", n
            ampo .+= eph_coefs[n], "Ntot(B1d+B1)", n
        end
        
        # phonon mode # 2
        if site_type == "Copper" && dim_copper_mode_2 > 1
            ampo .+= ω, "Nb2", n
            ampo .+= eph_coefs[n], "Ntot(B2d+B2)", n
        elseif site_type == "Oxygen" && dim_oxygen_mode_2 > 1
            ampo .+= ω, "Nb2", n
            ampo .+= eph_coefs[n], "Ntot(B2d+B2)", n
        end
        
        # phonon mode # 3
        if site_type == "Copper" && dim_copper_mode_3 > 1
            ampo .+= ω, "Nb3", n
            ampo .+= eph_coefs[n], "Ntot(B3d+B3)", n
        elseif site_type == "Oxygen" && dim_oxygen_mode_3 > 1
            ampo .+= ω, "Nb3", n
            ampo .+= eph_coefs[n], "Ntot(B3d+B3)", n
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
        ampo .-= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # oxygen-oxygen hopping 
    for b in pp_lattice
        ampo .-= b.sign*tpp, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign*tpp, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign*tpp, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign*tpp, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # electron-phonon nearest neighbour, oxygen-copper 
    for b in dp_lattice
        if site_labels[b.s1] == "Copper"
            # Put ph on copper, put e on oxygen 
            if dim_copper_mode_1 > 1
                ampo .+= g1dp, "B1dag+B1", b.s1, "Ntot", b.s2
            end
            if dim_copper_mode_2 > 1
                ampo .+= g1dp, "B2dag+B2", b.s1, "Ntot", b.s2
            end
            if dim_copper_mode_3 > 1
                ampo .+= g1dp, "B3dag+B3", b.s1, "Ntot", b.s2
            end
            # Put e on copper, put ph on oxygen 
            if dim_oxygen_mode_1 > 1
                ampo .+= g1pd, "Ntot", b.s1, "B1dag+B1", b.s2
            end
            if dim_oxygen_mode_2 > 1
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
        ampo .-= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .-= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .-= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .-= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1
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
                ampo .+= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2, "B1dag+B1", b.s2
                ampo .+= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1, "B1dag+B1", b.s2
                ampo .+= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B1dag+B1", b.s2
                ampo .+= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B1dag+B1", b.s2
            end

            # Mode 2 
            if dim_oxygen_mode_2 > 1
                ampo .+= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2, "B2dag+B2", b.s2
                ampo .+= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1, "B2dag+B2", b.s2
                ampo .+= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B2dag+B2", b.s2
                ampo .+= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B2dag+B2", b.s2
            end

            # Mode 3 
            if dim_oxygen_mode_3 > 1
                ampo .+= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2, "B3dag+B3", b.s2
                ampo .+= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1, "B3dag+B3", b.s2
                ampo .+= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B3dag+B3", b.s2
                ampo .+= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B3dag+B3", b.s2
            end

        elseif site_labels[b.s1] == "Oxygen"
            # Mode 1 (only add phonons to oxygen)
            if dim_oxygen_mode_1 > 1
                ampo .+= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2, "B1dag+B1", b.s1
                ampo .+= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1, "B1dag+B1", b.s1
                ampo .+= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B1dag+B1", b.s1
                ampo .+= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B1dag+B1", b.s1
            end

            # Mode 2 
            if dim_oxygen_mode_2 > 1
                ampo .+= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2, "B2dag+B2", b.s1
                ampo .+= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1, "B2dag+B2", b.s1
                ampo .+= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B2dag+B2", b.s1
                ampo .+= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B2dag+B2", b.s1
            end

            # Mode 3 
            if dim_oxygen_mode_3 > 1
                ampo .+= b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2, "B3dag+B3", b.s1
                ampo .+= b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1, "B3dag+B3", b.s1
                ampo .+= b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2, "B3dag+B3", b.s1
                ampo .+= b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1, "B3dag+B3", b.s1
            end
        end
    end

    return MPO(ampo,sites)
end

"""
Use this function for reconstructing from scratch given site indices 
"""
function ThreeBandModel(p::Parameters, sites::Vector{Index{Vector{Pair{QN, Int64}}}})
    H = make_ampo(p, sites)
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
    site_labels = make_coefficients(Nx+1, Ny, "Copper", "Oxygen", "Oxygen")
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
    expect(ψ,"Ntot")
end

function compute_entropy(ψ::MPS, b::Int)
    orthogonalize!(ψ, b)
    U,S,V = svd(ψ[b], (linkind(ψ, b-1), siteind(ψ,b)))
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
    return SvN
end

## DMRG ## 

function initialize_wavefcn(HM::ThreeBandModel, p::Parameters)
    # Extract the information 
    ups = "Up,"*string(p.init_phonons)
    downs = "Dn,"*string(p.init_phonons)
    emps = "Emp,"*string(p.init_phonons)
    
    state_arr = make_coefficients(p.Nx+1, p.Ny, ups, emps, emps)[1:p.Nsites]
    # Make every other hole a down 
    inds_arr = findall(x -> x==ups, state_arr)[begin:2:end]
    state_arr[inds_arr] .= downs
    # Last rung does not get any fermions
    state_arr[end-2*Ny+1:end] .= emps

    # Account for doping
    if p.doping > 0
        num_holes_doped = floor(Int, p.doping*p.N)
        @show num_holes_doped
        # put these holes on the px oxygen orbitals, spaced evenly apart 
        spacing = floor(Int,p.N/num_holes_doped) # spacing between unit cells 
        for (i,n) in enumerate(1:spacing:p.N)
            idx = to_site_number(p.Ny, n, "px")
            if isodd(i)
                state_arr[idx] = ups
            else
                state_arr[idx] = downs
            end
        end
    end

    num_holes_total = length(state_arr) - length(findall(x-> x==emps, state_arr))
    @show num_holes_total

    # NOTE: the QN of this state is preserved during DMRG
    productMPS(HM.sites,state_arr) 
end

function ladder_expectation(ϕ::MPS, opname::String, p::Parameters)
    density1D = expect(ϕ, opname)
    return reshape_into_lattice(density1D, p.Nx, p.Ny)
end

function expectation_tensor(T::ITensor, opname::String, s; normalization=1)
    val = scalar(T * op(opname, s) * dag(prime(T, s))) / normalization
end

function ladder_expectation_phonons(ϕ::MPS, p::Parameters)
    norm2_ϕ = norm(ϕ)^2
    sites = siteinds(ϕ)

    site_labels = make_coefficients(Nx+1, Ny, "Copper", "Oxygen", "Oxygen")
    density1D = zeros(p.Nsites, 3)
    for (n, site_type) in zip(1:p.Nsites, site_labels)
        # Phonon mode 1 
        if site_type == "Copper" && p.dim_copper_mode_1 > 1
            density1D[n,1] = expectation_tensor(ϕ[n], "Nb1", sites[n], normalization=norm2_ϕ)
        elseif site_type=="Oxygen" && p.dim_oxygen_mode_1 > 1
            density1D[n,1] = expectation_tensor(ϕ[n], "Nb1", sites[n], normalization=norm2_ϕ)
        end
        # Phonon mode 2
        if site_type == "Copper" && p.dim_copper_mode_2 > 1
            density1D[n,2] = expectation_tensor(ϕ[n], "Nb2", sites[n], normalization=norm2_ϕ)
        elseif site_type=="Oxygen" && p.dim_oxygen_mode_2 > 1
            density1D[n,2] = expectation_tensor(ϕ[n], "Nb2", sites[n], normalization=norm2_ϕ)
        end
        # Phonon mode 3 
        if site_type == "Copper" && p.dim_copper_mode_3 > 1
            density1D[n,3] = expectation_tensor(ϕ[n], "Nb3", sites[n], normalization=norm2_ϕ)
        elseif site_type=="Oxygen" && p.dim_oxygen_mode_3 > 1
            density1D[n,3] = expectation_tensor(ϕ[n], "Nb3", sites[n], normalization=norm2_ϕ)
        end
    end
    # reshape each dimension 
    density2D = zeros(p.Ny, p.Nx*3+2, 3)
    for i in 1:3
        density2D[:,:,i] = reshape_into_lattice(density1D[:,i], p.Nx, p.Ny)
    end 
    return density2D 
end

function get_sweeps(p::Parameters, DMRG_numsweeps_per_save::Int)
    # Set sweep params
    sweeps = Sweeps(p.DMRG_numsweeps)
    setnoise!(sweeps, p.DMRG_noise...) # Very important to use noise for this model
    setmaxdim!(sweeps, p.DMRG_maxdim...)
    setcutoff!(sweeps, p.DMRG_cutoff...)
    nsweep = p.DMRG_numsweeps
    noise = sweeps.noise
    maxdim = sweeps.maxdim
    cutoff = sweeps.cutoff 

    if !isnothing(DMRG_numsweeps_per_save)
        sweeps = Sweeps(DMRG_numsweeps_per_save)
        setnoise!(sweeps, noise[1:DMRG_numsweeps_per_save]...)
        setmaxdim!(sweeps, maxdim[1:DMRG_numsweeps_per_save]...)
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
        setcutoff!(sweeps, cutoff...)
        nsweep += DMRG_numsweeps_per_save
    end

    return sweeps, nsweep, noise, maxdim, cutoff 
end

function run_DMRG(dmrg_results::DMRGResults, HM::ThreeBandModel, p::Parameters; 
                    DMRG_numsweeps_per_save::Union{Nothing,Int}=nothing, 
                    alg="divide_and_conquer", disk_save=false)
    # Set DMRG params
    sweeps, nsweep, noise, maxdim, cutoff = get_sweeps(dmrg_results, DMRG_numsweeps_per_save) 
    
    # Load in the last wavefunction 
    ϕ0 = dmrg_results.ground_state
    if p.DMRG_LBO # If performing local basis optimization
        @assert 1==0 "Need to load in the old Rs from a previous run"
        energy, ϕ, Rs = dmrg_lbo(HM.mpo, ϕ0, sweeps, alg=alg, LBO=true, 
                                    max_LBO_dim=p.max_LBO_dim, min_LBO_dim=p.min_LBO_dim)
    else
        if disk_save
            energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg, write_when_maxdim_exceeds=p.DMRG_maxdim)
        else
            energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg)
        end
        Rs = 0
    end

    # compute some quantities
    entropy = compute_entropy(ϕ, p.mid)
    charge_density = ladder_expectation(ϕ, "Ntot", p)
    spin_density = ladder_expectation(ϕ, "Sz", p)
    phonon_density = ladder_expectation_phonons(ϕ, p)
    return DMRGResults(nsweep, maxdim, cutoff, noise, ϕ, energy, entropy, 
                        Rs, charge_density, phonon_density, spin_density)
end

function run_DMRG(HM::ThreeBandModel, p::Parameters; 
                    DMRG_numsweeps_per_save::Union{Nothing,Int}=nothing, 
                    alg="divide_and_conquer", disk_save=false)
    
    # Set DMRG params
    sweeps, nsweep, noise, maxdim, cutoff = get_sweeps(p, DMRG_numsweeps_per_save)

    # Initialize the wavefunction 
    ϕ0 = initialize_wavefcn(HM,p)
    
    # Run DMRG 
    if p.DMRG_LBO # If performing local basis optimization
        energy, ϕ, Rs = dmrg_lbo(HM.mpo, ϕ0, sweeps, alg=alg, LBO=true, 
                                    max_LBO_dim=p.max_LBO_dim, min_LBO_dim=p.min_LBO_dim)
    else
        if disk_save
            energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg, write_when_maxdim_exceeds=p.DMRG_maxdim)
        else
            energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg)
        end
        Rs = 0
    end

    # compute some quantities
    entropy = compute_entropy(ϕ, p.mid)
    charge_density = ladder_expectation(ϕ, "Ntot", p)
    spin_density = ladder_expectation(ϕ, "Sz", p)
    phonon_density = ladder_expectation_phonons(ϕ, p)

    return DMRGResults(nsweep, maxdim, cutoff, noise, ϕ, energy, entropy, 
                        Rs, charge_density, phonon_density, spin_density)
end

## CORRELATION FUNCTIONS ## 

function apply_onesite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
    ϕ = copy(ϕ) 

    ## Account for fermion sign using Jordan-Wigner strings ##
    if opname == "Cup" || opname == "Cdn"
        ϕ = apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
        for i in reverse(1:(siteidx-1)) # Don't act with string on-site
            ϕ = apply_op(ϕ, op("F",sites[i]), i)
        end
        return ϕ

    elseif opname == "Cdagup" || opname == "Cdagdn"
        for i in 1:(siteidx-1) # Don't act with string on-site
            ϕ = apply_op(ϕ, op("F",sites[i]), i)
        end
        ϕ = apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
        return ϕ
    end

    # Otherwise, just apply the operator as usual
    return apply_op(ϕ, op(opname,sites[siteidx]), siteidx)
end

function Π̂(sites, s1::Int, s2::Int)
    # Π = 1/sqrt(2)*(op("Cup",sites[s1])*op("Cdn",sites[s2])
    #                     + op("Cdn",sites[s1])*op("Cup",sites[s2]))

    if s1 > s2 # ordering matters! these are fermions we're talking about
        @warn "s1 > s2"
        s1, s2 = s2, s1
    end

    # First part 
    Π1 = op("Cup",sites[s1])
    for s in (s1+1):(s2-1)
        Π1 *= op("F",sites[s])
    end
    Π1 *= op("Cdn",sites[s2])

    # Second part 
    Π2 = op("Cdn",sites[s1])
    for s in (s1+1):(s2-1)
        Π2 *= op("F",sites[s])
    end
    Π2 *= op("Cup",sites[s2])

    # Take sum and normalize
    Π = 1/sqrt(2)*(Π1 + Π2)
end

function Δ̂(sites, s1::Int, s2::Int)
    # Δ = 1/sqrt(2)*(op("Cup",sites[s1])*op("Cdn",sites[s2])
    #                         - op("Cdn",sites[s1])*op("Cup",sites[s2]))

    if s1 > s2 # ordering matters! these are fermions we're talking about
        @warn "s1 > s2"
        s1, s2 = s2, s1
    end

    # First part
    Δ1 = op("Cup",sites[s1])
    for s in (s1+1):(s2-1)
        Δ1 *= op("F",sites[s])
    end
    Δ1 *= op("Cdn",sites[s2])

    # Second part 
    Δ2 = op("Cdn",sites[s1])
    for s in (s1+1):(s2-1)
        Δ2 *= op("F",sites[s])
    end
    Δ2 *= op("Cup",sites[s2])

    # Take difference and normalize
    Δ = 1/sqrt(2)*(Δ1 - Δ2)
end

function apply_twosite_operator(ϕ::MPS, opname::String, sites, s1::Int, s2::Int)
    ϕ = copy(ϕ)
    if opname == "pSC"
        Π = Π̂(sites, s1, s2)
        return apply(Π, ϕ)     
    elseif opname == "dSC"
        Δ = Δ̂(sites, s1, s2)
        return apply(Δ, ϕ)
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
                                            buffer=1)

    # Compute spin and charge correlations
    println("Computing spin correlations for dx-dx bond")
    start, stop, spin_corr = compute_equilibrium_onsite_correlation(dmrg_results,HM,p,"dx-dx","spin")

    println("Computing charge correlations for dx-dx bond")
    _,_,charge_corr = compute_equilibrium_onsite_correlation(dmrg_results,HM,p,"dx-dx","charge")

    # Compute d-wave pairfield correlations 
    pairfield_corrs = compute_equilibrium_pairfield_correlations(dmrg_results, HM, p)
    
    return EquilibriumCorrelations(start, stop, spin_corr, charge_corr, pairfield_corrs...)
end

function compute_equilibrium_pairfield_correlations(dmrg_results::DMRGResults, 
                                                HM::ThreeBandModel, p::Parameters)

    # There are 6 pair-field correlations to compute 
    println("Computing dSC correlations for dx-dx bond")
    _,_,dxdx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dx-dx", "dx-dx", "dSC")
    println("Computing dSC correlations for d-px bond")
    _,_,dpx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "d-px", "d-px", "dSC")
    println("Computing dSC correlations for dy-dy bond")
    _,_,dydy = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dy-dy", "dy-dy", "dSC")
    println("Computing dSC correlations for py-d bond")
    _,_,pyd = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-d", "py-d", "dSC")
    println("Computing dSC correlations for py-px bond")
    _,_,pypx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-px", "py-px", "dSC")
    println("Computing dSC correlations for py1-px2 bond")
    _,_,py1px2 = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py1-px2", "py1-px2", "dSC")

    return [dxdx, dpx, dydy, pyd, pypx, py1px2] # transpose for compatibiity w plotting fcn
end

function get_bonds(bondtype::String, lattice_indices; row=1)
    if bondtype=="d-px" 
        s1 = 2
        s2 = 3
        r1, r2 = row, row
    elseif bondtype=="px-d"
        s1 = 3
        s2 = 2
        r1, r2 = row, row
    elseif bondtype=="dx-dx"
        s1 = 2
        s2 = 5
        r1, r2 = row, row
    elseif bondtype=="dy-dy"
        s1 = 2
        s2 = 2
        r1, r2 = 1, 2
    elseif bondtype=="py-d"
        s1 = 1
        s2 = 2
        r1, r2 = row, row
    elseif bondtype=="py1-px2"
        s1 = 1
        s2 = 3
        r1, r2 = 1, 2
    elseif bondtype=="py-px"
        s1 = 1
        s2 = 3
        r1, r2 = row, row
    else
        @error "Bond type not recognized"
        return nothing
    end

    refbond = (lattice_indices[r1, s1],lattice_indices[r2, s2])
    sites_1 = lattice_indices[r1, s1:3:end] 
    sites_2 = lattice_indices[r2, s2:3:end] 
    bonds = [(sites_1[n], sites_2[n]) for n in 1:length(sites_2)]

    return refbond, bonds
end

function compute_equilibrium_onsite_correlation(dmrg_results::DMRGResults, 
                                                HM::ThreeBandModel, p::Parameters,
                                                bondtype::String, corrtype::String;
                                                buffer=1, row=1)
    Nx, Ny = p.Nx, p.Ny
    ϕ = copy(dmrg_results.ground_state)
    Nsites = length(ϕ)
    start = 4 + (buffer-1)*3 # discard buffer # of unit cells 
    stop = 3*Nx - (buffer-1)*3 - 2 # discard the last rung and buffer # of unit cells

    # Compute the correlations for each row separately 
    lattice_indices = reshape_into_lattice(collect(1:Nsites), Nx, Ny)
    lattice_indices = convert.(Int, lattice_indices) 
    lattice_indices = lattice_indices[:,start:stop] 

    # Compute the equilibrium correlations within this chain 
    _, bonds = get_bonds(bondtype, lattice_indices, row=row)
    corr = onsite_correlation(ϕ, bonds, corrtype, HM.sites)

    return start,stop,corr 
end

function compute_equilibrium_pairfield_correlation(dmrg_results::DMRGResults, 
                                                HM::ThreeBandModel, p::Parameters,
                                                bond1::String, bond2::String, SCtype::String;
                                                buffer=1, row=1)
    Nx, Ny = p.Nx, p.Ny
    ϕ = copy(dmrg_results.ground_state)
    Nsites = length(ϕ)
    start = 4 + (buffer-1)*3 # discard buffer # of unit cells 
    stop = 3*Nx - (buffer-1)*3 - 2 # discard the last rung and buffer # of unit cells

    # Compute the correlations for each row separately 
    lattice_indices = reshape_into_lattice(collect(1:Nsites), Nx, Ny)
    lattice_indices = convert.(Int, lattice_indices) 
    lattice_indices = lattice_indices[:,start:stop] 

    # Compute the equilibrium correlations within this chain 
    refbond, _ = get_bonds(bond1, lattice_indices, row=row)
    _, bonds = get_bonds(bond2, lattice_indices, row=row)
    corr = bond_correlation(ϕ, refbond, bonds, SCtype, HM.sites)

    return start,stop,corr 
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

    if corrtype=="spin"
        ψ = apply_onesite_operator(ϕ, "Sz", sites, j)
        function compute_corr_spin(i::Int)
            Szψ = apply_onesite_operator(ψ, "Sz", sites, i)
            return inner(ϕ,Szψ)
        end
        return compute_corr_spin.(indices)
    elseif corrtype=="charge"
        ψ = apply_onesite_operator(ϕ, "Ntot", sites, j)
        function compute_corr_charge(i::Int)
            Ntotψ = apply_onesite_operator(ψ, "Ntot", sites, i)
            return inner(ϕ,Ntotψ)
        end
        ninj = compute_corr_charge.(indices)
        ni = expect(ϕ, "Ntot")
        nj = ni[j]
        return ninj - nj .* (ni[indices])
    elseif corrtype=="sSC"
        ψ = apply_onesite_operator(ϕ, "Cupdn", sites, j)
        function compute_corr_sSC(i::Int)
            Σ_iϕ = apply_onesite_operator(ϕ, "Cupdn", sites, i)
            return inner(ψ,Σ_iϕ)
        end
        return compute_corr_sSC.(indices)
    end
    @error "Unknown correlation type"
end

function bond_correlation(ϕ::MPS, refbond, bonds, corrtype::String, sites)
    if corrtype=="sSC"
        indices = unzip(bonds) # we don't actually care about bonds, bc it's all on single sites
        ψ = apply_onesite_operator(ϕ, "Cupdn", sites, indices[1])
        function compute_corr_sSC(i::Int)
            Σ_iϕ = apply_onesite_operator(ϕ, "Cupdn", sites, i)
            return inner(ψ,Σ_iϕ)
        end
        return compute_corr_sSC.(indices)

    elseif corrtype=="pSC"
        ψ = apply_twosite_operator(ϕ, "pSC", sites, refbond[1], refbond[2])
        function compute_corr_pSC(bond)
            Π_iϕ = apply_twosite_operator(ϕ, "pSC", sites, bond[1], bond[2])
            return inner(ψ,Π_iϕ)
        end
        return compute_corr_pSC.(bonds)

    elseif corrtype=="dSC"
        ψ = apply_twosite_operator(ϕ, "dSC", sites, refbond[1], refbond[2])
        function compute_corr_dSC(bond)
            Δ_iϕ = apply_twosite_operator(ϕ, "dSC", sites, bond[1], bond[2])
            return inner(ψ,Δ_iϕ)
        end
        return compute_corr_dSC.(bonds)        
    end
    @error "Unknown correlation type"
end