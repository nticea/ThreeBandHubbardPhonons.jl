using ITensors
include("cuprate_lattice.jl")

struct Parameters
    # Model 
    N::Int # number of unit cells
    Nx::Int # number of unit cells along x
    Ny::Int # number of unti cells along y 
    Nsites::Int # number of sites in the MPS representation 
    yperiodic::Bool
    doping::Real
    max_phonons::Int
    init_phonons::Int

    # Hamiltonian
    μ::Real
    εd::Real
    εp::Real
    tpd::Real
    tpp::Real
    Upd::Real
    Upp::Real
    Udd::Real
    ω::Real
    g0pp::Real
    g0dd::Real
    g1pd::Real
    g1dp::Real
    g1pp::Real

    # DMRG parameters
    DMRG_numsweeps
    DMRG_noise
    DMRG_maxdim
    DMRG_cutoff

    # LBO parameters 
    DMRG_LBO::Bool
    max_LBO_dim::Int
    min_LBO_dim::Int

    # TEBD parameters
    mid::Int
    T::Int
    τ::Real
    TEBD_cutoff
    TEBD_maxdim
    TEBD_LBO::Bool
end

struct ThreeBandModel
    # Lattice
    sites 
    # Hamiltonian (in two forms)
    mpo::MPO
    gates 
end

struct DMRGResults
    ground_state
    ground_state_energy
    ground_state_entropy 
    optimized_basis
    charge_density
    phonon_density
    spin_density
end

struct EquilibriumCorrelations
    start
    stop
    spin
    charge
    sSC
    pSC
    dSC
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
    max_phonons::Int=1, 
    init_phonons::Int=0,

    μ=0,
    εd=0,
    εp=3,
    tpd=1,
    tpp=0.3,
    Upd=0.5,
    Upp=3,
    Udd=8,
    ω::Real,
    g0pp::Real,
    g0dd::Real,
    g1pd::Real,
    g1dp::Real,
    g1pp::Real,

    DMRG_numsweeps::Int=20, DMRG_noise=nothing, 
    DMRG_maxdim=nothing, DMRG_cutoff=nothing, DMRG_LBO=false,
    max_LBO_dim=nothing, min_LBO_dim=4,
    T::Int=25, τ::Real=0.1, TEBD_cutoff=1E-14, TEBD_maxdim=400, TEBD_LBO=false)

    if isnothing(ω)
        ω = 0.5*t
    end
    if isnothing(g0pp)
        g0pp = 0.05*t
    end
    if isnothing(g0dd)
        g0dd = 0.05*t
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

    if max_phonons==0
        @assert ω==g0==g1==0
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
    Nsites = N*3 # 3 sites per unit cell 
    mid = ceil(Int,Nsites/2) # midpoint of the DMRG chain 

    return Parameters(N, Nx, Ny, Nsites, yperiodic, doping, max_phonons, init_phonons,
            μ, εd, εp, tpd, tpp, Upd, Upp, Udd, ω, g0pp, g0dd, g1pd, g1dp, g1pp,
            DMRG_numsweeps, DMRG_noise, DMRG_maxdim, DMRG_cutoff, DMRG_LBO,
            max_LBO_dim, min_LBO_dim, mid, T, τ, TEBD_cutoff, TEBD_maxdim, TEBD_LBO)
    
end

function make_ampo(p::Parameters, sites::Vector{Index{Vector{Pair{QN, Int64}}}})
    N, Nsites, Nx, Ny, max_phonons, yperiodic = p.N, p.Nsites, p.Nx, p.Ny, p.max_phonons, p.yperiodic
    μ, εd, εp, tpd, tpp, Upd, Upp, Udd, ω, g0pp, g0dd, g1pd, g1dp, g1pp = p.μ, p.εd, p.εp, p.tpd, p.tpp, p.Upd, p.Upp, p.Udd, p.ω, p.g0pp, p.g0dd, p.g1pd, p.g1dp, p.g1pp
    dp_lattice = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic)
    pp_lattice = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic)

    # make the hamiltonian 
    ampo = OpSum()

    # bond terms
    for b in lattice
        ampo .+= -t, "Cdagup", b.s1, "Cup", b.s2
        ampo .+= -t, "Cdagup", b.s2, "Cup", b.s1
        ampo .+= -t, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .+= -t, "Cdagdn", b.s2, "Cdn", b.s1
    end

    # on-site terms 
    for j=1:N-1       
        # ∑_j U * n_j↑ n_j↓ 
        ampo += U,"Nupdn",j,"I",j

        # ∑_j ω * nb_j
        if max_phonons >= 1
            ampo += ω,"Nb",j

            # # ∑_j g0 * nf_j (b^†_j + b_j)
            ampo += g0,"Ntot(Bd+B)",j

            # # ∑_⟨ij⟩ g1 * nf_j (b^†_i + b_i)
            ampo += g1,"Ntot",j,"Bdag+B",j+1
            ampo += g1,"Ntot",j+1,"Bdag+B",j

            # quartic term
            ampo += λ,"Nb^2",j
        end
    end
    # Edge site
    ampo += U,"Nupdn",N
    if max_phonons >= 1
        ampo += ω,"Nb",N
        ampo += g0,"Ntot(Bd+B)",N
        ampo += λ,"Nb^2",N
    end
    return MPO(ampo,sites)
end