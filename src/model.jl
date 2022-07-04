using ITensors
using Statistics: mean
include("cuprate_lattice.jl")
include("site_hubbholst.jl")

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

    if max_phonons==0
        @assert ω==g0pp==g0dd==g1pd==g1dp==g1pp==0
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
    Nsites, Nx, Ny, yperiodic, max_phonons = p.Nsites, p.Nx, p.Ny, p.yperiodic, p.max_phonons
    μ, εd, εp, tpd, tpp, Upd, Upp, Udd, ω, g0pp, g0dd, g1pd, g1dp, g1pp = p.μ, p.εd, p.εp, p.tpd, p.tpp, p.Upd, p.Upp, p.Udd, p.ω, p.g0pp, p.g0dd, p.g1pd, p.g1dp, p.g1pp
    dp_lattice = OxygenCopper_lattice(Nx, Ny; yperiodic=yperiodic)
    pp_lattice = OxygenOxygen_lattice(Nx, Ny; yperiodic=yperiodic)

    # make the hamiltonian 
    ampo = OpSum()

    # on-site terms 
    μ_coefs = make_coefficients(Nx, Ny, μ-εd, μ-εp, μ-εp)
    U_coefs = make_coefficients(Nx, Ny, Udd, Upp, Upp)
    eph_coefs = make_coefficients(Nx, Ny, g0dd, g0pp, g0pp)
    for n in 1:Nsites
        # chemical potential term
        ampo .+= μ_coefs[n], "Ntot", n
        # on-site repulsion
        ampo .+= U_coefs[n], "Nup", n, "Ndn", n
        if max_phonons>0
            # phonon 
            ampo .+= ω, "Nb", n
            # on-site e-ph interactions
            ampo .+= eph_coefs[n], "Ntot(Bd+B)", n
        end
    end

    # repulsion copper-oxygen
    for b in dp_lattice
        ampo .+= Upd, "Nup", b.s1, "Nup", b.s2
        ampo .+= Upd, "Ndn", b.s1, "Ndn", b.s2
        ampo .+= Upd, "Nup", b.s1, "Ndn", b.s2
        ampo .+= Upd, "Ndn", b.s1, "Nup", b.s2
    end

    # copper-oxygen hopping 
    for b in dp_lattice
        ampo .+= -b.sign*tpd, "Cdagup", b.s1, "Cup", b.s2
        ampo .+= -b.sign*tpd, "Cdagup", b.s2, "Cup", b.s1
        ampo .+= -b.sign*tpd, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .+= -b.sign*tpd, "Cdagdn", b.s2, "Cdn", b.s1
    end
    # oxygen-oxygen hopping 
    for b in pp_lattice
        ampo .+= -b.sign*tpp, "Cdagup", b.s1, "Cup", b.s2
        ampo .+= -b.sign*tpp, "Cdagup", b.s2, "Cup", b.s1
        ampo .+= -b.sign*tpp, "Cdagdn", b.s1, "Cdn", b.s2
        ampo .+= -b.sign*tpp, "Cdagdn", b.s2, "Cdn", b.s1
    end

    if max_phonons>0
        # electron-phonon nearest neighbour
        for b in dp_lattice
            # phonon on copper, electron on oxygen
            ampo .+= g1dp, "Bdag+B", b.s1, "Ntot", b.s2 
            # phonon on oxygen, electron on copper
            ampo .+= g1pd, "Bdag+B", b.s2, "Ntot", b.s1 
        end
        for b in pp_lattice
            # phonon on oxygen, electron on oxygen
            ampo .+= g1pp, "Bdag+B", b.s1, "Ntot", b.s2
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

"""
    ThreeBandModel(p::Parameters)
Make a three-band Hubbard model with phonons given a set of input parameters 
"""
function ThreeBandModel(p::Parameters)
    sites = siteinds("HubHolst", p.Nsites; dim=p.max_phonons+1)
    return ThreeBandModel(p, sites)
end

## GENERIC STATE FUNCTIONS ##

function compute_overlap(ψ1::MPS, ψ2::MPS)
    LinearAlgebra.norm(inner(ψ1, ψ2))
end

function compute_phonon_number(ψ::MPS)
    expect(ψ,"Nb")
end

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
    
    state = [isodd(n) ? ups : downs for n=1:p.Nsites] 

    # Account for doping
    if p.doping > 0
        spacing = floor(Int,1/p.doping)
        state[1:spacing:end] .= emps
    end

    # NOTE: the QN of this state is preserved during DMRG
    productMPS(HM.sites,state) 
end

function ladder_expectation(ϕ::MPS, opname::String, p::Parameters)
    density1D = expect(ϕ, opname)
    return reshape_into_lattice(density1D, p.Nx, p.Ny)
end

function run_DMRG(HM::ThreeBandModel, p::Parameters; alg="divide_and_conquer")
    # Set DMRG params
    sweeps = Sweeps(p.DMRG_numsweeps)
    setnoise!(sweeps, p.DMRG_noise...) # Very important to use noise for this model
    setmaxdim!(sweeps, p.DMRG_maxdim...)
    setcutoff!(sweeps, p.DMRG_cutoff...) 
    
    ϕ0 = initialize_wavefcn(HM,p)
    @show flux(ϕ0)
    if p.DMRG_LBO # If performing local basis optimization
        @error "LBO not implemented yet"
        # energy, ϕ, Rs = dmrg_lbo(HM.mpo, ϕ0, sweeps, alg=alg, LBO=true, 
        #                             max_LBO_dim=p.max_LBO_dim, min_LBO_dim=p.min_LBO_dim)
    else
        energy, ϕ = dmrg(HM.mpo, ϕ0, sweeps, alg=alg)
        Rs = nothing
    end

    # compute some quantities
    entropy = compute_entropy(ϕ, p.mid)
    charge_density = ladder_expectation(ϕ, "Ntot", p)
    spin_density = ladder_expectation(ϕ, "Sz", p)
    if p.max_phonons > 0
        phonon_density = ladder_expectation(ϕ, "Nb", p)
    else
        phonon_density = zeros(p.Nx, p.Ny)
    end

    return DMRGResults(ϕ, energy, entropy, Rs, charge_density, 
                            phonon_density, spin_density)
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

function apply_twosite_operator(ϕ::MPS, opname::String, sites, siteidx::Int)
    ϕ = copy(ϕ)
    if opname == "pSC"
        # Apply F string first 
        ϕ = apply_op(ϕ, op("F",sites[siteidx]), siteidx)

        # Then apply the operator
        Π = 1/sqrt(2)*(op("Cup",sites[siteidx])*op("Cdn",sites[siteidx+1])
                            + op("Cdn",sites[siteidx])*op("Cup",sites[siteidx+1]))
        return apply_op_twosite(ϕ,Π,siteidx)
                    
    elseif opname == "pSCdag"
        # Apply the operator first
        Π_dag = 1/sqrt(2)*(op("Cdagdn",sites[siteidx+1])*op("Cdagup",sites[siteidx])
                            + op("Cdagup",sites[siteidx+1])*op("Cdagdn",sites[siteidx]))
        ϕ = apply_op_twosite(ϕ,Π_dag,siteidx)

        # Apply the F string 
        return apply_op(ϕ, op("F",sites[siteidx]), siteidx)

    elseif opname == "dSC"
        # Apply F string first 
        ϕ = apply_op(ϕ, op("F",sites[siteidx]), siteidx)

        # Then apply the operator
        Δ = 1/sqrt(2)*(op("Cup",sites[siteidx])*op("Cdn",sites[siteidx+1])
                            - op("Cdn",sites[siteidx])*op("Cup",sites[siteidx+1]))
        return apply_op_twosite(ϕ,Δ,siteidx)

    elseif opname == "dSCdag"
        # Apply the operator first
        Δ_dag = 1/sqrt(2)*(op("Cdagdn",sites[siteidx+1])*op("Cdagup",sites[siteidx])
                            - op("Cdagup",sites[siteidx+1])*op("Cdagdn",sites[siteidx]))
        ϕ = apply_op_twosite(ϕ,Δ_dag,siteidx)

        # Apply the F string 
        return apply_op(ϕ, op("F",sites[siteidx]), siteidx)
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

function apply_op_twosite(ϕ::MPS, G::ITensor, siteidx::Int; cutoff=1E-8)
    ϕ = copy(ϕ)
    # Note: siteidx corresponds to the leftermost site
    orthogonalize!(ϕ,siteidx)
    wf = (ϕ[siteidx] * ϕ[siteidx+1]) * G
    noprime!(wf)

    inds_site = uniqueinds(ϕ[siteidx],ϕ[siteidx+1])
    U,S,V = svd(wf,inds_site,cutoff=cutoff)
    ϕ[siteidx] = U
    ϕ[siteidx+1] = S*V
    return ϕ
end

function compute_all_equilibrium_correlations(dmrg_results::DMRGResults, 
                                            HM::ThreeBandModel, p::Parameters;
                                            start=nothing, stop=nothing)
    ϕ = copy(dmrg_results.ground_state)

    # Avoid edge effects by only considering the middle 
    # N = p.Nsites
    # if isnothing(start)
    #     start = floor(Int,0.25*N)
    # end 
    # if isnothing(stop)
    #     stop = ceil(Int,0.75*N)
    # end
    # lencorrs = stop-start+1
    # if (stop-start+1)%Ny != 0
    #     stop += (Ny - lencorrs%Ny)
    # end

    # @assert stop <= N
    # @assert (stop-start+1)%Ny == 0

    # Compute the correlations for each row separately 
    lattice_indices = reshape_into_lattice(collect(1:p.Nsites), p.Nx, p.Ny)
    lattice_indices = convert.(Int, lattice_indices)

    # Iterate through all the correlations types and compute them 
    corrtypes = ["spin","charge","sSC","pSC","dSC"]
    corrs = []
    for corrtype in corrtypes
        corr = zeros(Ny, 3*Nx)
        # discard the first unit cell due to edge effects 
        for y in 1:Ny
            if Ny ==1
                indices = lattice_indices[y,4:end-2]
                corr[y,4:end-2] = equilibrium_correlations(ϕ,indices,corrtype,HM.sites)
            elseif y == 1
                indices = lattice_indices[y,4:end]
                corr[y,4:end] = equilibrium_correlations(ϕ,indices,corrtype,HM.sites)
            elseif y == Ny
                indices = lattice_indices[y,1:end-2]
                corr[y,1:end-2] = equilibrium_correlations(ϕ,indices,corrtype,HM.sites)
            else
                indices = lattice_indices[y,:]
                corr[y,:] = equilibrium_correlations(ϕ,indices,corrtype,HM.sites)
            end
        end
        # Fill in the first 3 and last 3 elements by taking averages from the other ladders
        if Ny > 1
            corr[1,1:3] = mean(corr[2:end,1:3], dims=1)
            corr[Ny,end-2:end] = mean(corr[1:Ny-1,end-2:end], dims=1)
        else # we trim 
            corr = corr[:, 4:end-2]
        end

        push!(corrs,corr)
    end

    EquilibriumCorrelations(start, stop, corrs...)
end

function equilibrium_correlations(ϕ::MPS, indices, corrtype::String,sites)
    ϕ = copy(ϕ)
    j = indices[1]
    if corrtype=="spin"
        return correlation_matrix(ϕ, "Sz", "Sz")[indices,j]
    elseif corrtype=="charge"
        ninj = correlation_matrix(ϕ, "Ntot", "Ntot")[indices,j]
        ni = expect(ϕ, "Ntot")
        nj = ni[j]
        return ninj - nj .* (ni[indices])
    elseif corrtype=="sSC"
        ψ = apply_onesite_operator(ϕ, "Cupdn", sites, j)
        function compute_corr_sSC(i::Int)
            Σ_iψ = apply_onesite_operator(ψ, "Cdagupdn", sites, i)
            return inner(ϕ,Σ_iψ)
        end
        return compute_corr_sSC.(indices)
    elseif corrtype=="pSC"
        ψ = apply_twosite_operator(ϕ, "pSC", sites, j)
        function compute_corr_pSC(i::Int)
            Π_iψ = apply_twosite_operator(ψ, "pSCdag", sites, i)
            return inner(ϕ,Π_iψ)
        end
        return compute_corr_pSC.(indices)
    elseif corrtype=="dSC"
        ψ = apply_twosite_operator(ϕ, "dSC", sites, j)
        function compute_corr_dSC(i::Int)
            Δ_iψ = apply_twosite_operator(ψ, "dSCdag", sites, i)
            return inner(ϕ,Δ_iψ)
        end
        return compute_corr_dSC.(indices)        
    end
    @error "Unknown correlation type"
end