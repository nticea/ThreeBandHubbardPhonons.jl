using ITensors
using ITensors: AbstractMPS
using TimerOutputs
using Printf
using KrylovKit: eigsolve

function permute(
  M::AbstractMPS, ::Tuple{typeof(linkind),typeof(siteinds),typeof(linkind)}
)::typeof(M)
  M̃ = typeof(M)(length(M))
  for n in 1:length(M)
    lₙ₋₁ = linkind(M, n - 1)
    lₙ = linkind(M, n)
    s⃗ₙ = sort(Tuple(siteinds(M, n)); by=plev)
    M̃[n] = ITensors.permute(M[n], filter(!isnothing, (lₙ₋₁, s⃗ₙ..., lₙ)))
  end
  set_ortho_lims!(M̃, ortho_lims(M))
  return M̃
end

"""
    dmrg(H::MPO,psi0::MPS,sweeps::Sweeps;kwargs...)
                    
Use the density matrix renormalization group (DMRG) algorithm
to optimize a matrix product state (MPS) such that it is the
eigenvector of lowest eigenvalue of a Hermitian matrix H,
represented as a matrix product operator (MPO).
The MPS `psi0` is used to initialize the MPS to be optimized,
and the `sweeps` object determines the parameters used to 
control the DMRG algorithm.

Returns:
* `energy::Float64` - eigenvalue of the optimized MPS
* `psi::MPS` - optimized MPS

Optional keyword arguments:
* `outputlevel::Int = 1` - larger outputlevel values make DMRG print more information and 0 means no output
* `observer` - object implementing the [Observer](@ref observer) interface which can perform measurements and stop DMRG early
* `write_when_maxdim_exceeds::Int` - when the allowed maxdim exceeds this value, begin saving tensors to disk to free memory in large calculations
"""
function dmrg_lbo(H::MPO, psi0::MPS, sweeps::Sweeps; kwargs...)
  ITensors.check_hascommoninds(siteinds, H, psi0)
  ITensors.check_hascommoninds(siteinds, H, psi0')
  # Permute the indices to have a better memory layout
  # and minimize permutations
  H = permute(H, (linkind, siteinds, linkind))
  PH = ProjMPO(H)
  return dmrg_lbo(PH, psi0, sweeps; kwargs...)
end

"""
    dmrg(Hs::Vector{MPO},psi0::MPS,sweeps::Sweeps;kwargs...)
                    
Use the density matrix renormalization group (DMRG) algorithm
to optimize a matrix product state (MPS) such that it is the
eigenvector of lowest eigenvalue of a Hermitian matrix H.
The MPS `psi0` is used to initialize the MPS to be optimized,
and the `sweeps` object determines the parameters used to 
control the DMRG algorithm.

This version of `dmrg` accepts a representation of H as a
Vector of MPOs, Hs = [H1,H2,H3,...] such that H is defined
as H = H1+H2+H3+...
Note that this sum of MPOs is not actually computed; rather
the set of MPOs [H1,H2,H3,..] is efficiently looped over at 
each step of the DMRG algorithm when optimizing the MPS.

Returns:
* `energy::Float64` - eigenvalue of the optimized MPS
* `psi::MPS` - optimized MPS
"""
function dmrg_lbo(Hs::Vector{MPO}, psi0::MPS, sweeps::Sweeps; kwargs...)
  for H in Hs
    ITensors.check_hascommoninds(siteinds, H, psi0)
    ITensors.check_hascommoninds(siteinds, H, psi0')
  end
  Hs .= permute.(Hs, Ref((linkind, siteinds, linkind)))
  PHS = ProjMPOSum(Hs)
  return dmrg_lbo(PHS, psi0, sweeps; kwargs...)
end

"""
    dmrg(H::MPO,Ms::Vector{MPS},psi0::MPS,sweeps::Sweeps;kwargs...)
                    
Use the density matrix renormalization group (DMRG) algorithm
to optimize a matrix product state (MPS) such that it is the
eigenvector of lowest eigenvalue of a Hermitian matrix H,
subject to the constraint that the MPS is orthogonal to each
of the MPS provided in the Vector `Ms`. The orthogonality
constraint is approximately enforced by adding to H terms of 
the form w|M1><M1| + w|M2><M2| + ... where Ms=[M1,M2,...] and
w is the "weight" parameter, which can be adjusted through the
optional `weight` keyword argument.
The MPS `psi0` is used to initialize the MPS to be optimized,
and the `sweeps` object determines the parameters used to 
control the DMRG algorithm.

Returns:
* `energy::Float64` - eigenvalue of the optimized MPS
* `psi::MPS` - optimized MPS
"""
function dmrg_lbo(
  H::MPO, Ms::Vector{MPS}, psi0::MPS, sweeps::Sweeps; kwargs...)
  ITensors.check_hascommoninds(siteinds, H, psi0)
  ITensors.check_hascommoninds(siteinds, H, psi0')
  for M in Ms
    ITensors.check_hascommoninds(siteinds, M, psi0)
  end
  H = permute(H, (linkind, siteinds, linkind))
  Ms .= permute.(Ms, Ref((linkind, siteinds, linkind)))
  weight = get(kwargs, :weight, 1.0)
  PMM = ProjMPO_MPS(H, Ms; weight=weight)
  return dmrg_lbo(PMM, psi0, sweeps; kwargs...)
end

function dmrg_lbo(PH, psi0::MPS, sweeps::Sweeps; kwargs...)
    if length(psi0) == 1
      error(
        "`dmrg` currently does not support system sizes of 1. You can diagonalize the MPO tensor directly with tools like `LinearAlgebra.eigen`, `KrylovKit.eigsolve`, etc.",
      )
    end
  
    ITensors.@debug_check begin
      # Debug level checks
      # Enable with ITensors.enable_debug_checks()
      checkflux(psi0)
      checkflux(PH)
    end
  
    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    obs = get(kwargs, :observer, NoObserver())
    outputlevel::Int = get(kwargs, :outputlevel, 1)
    LBO::Bool = get(kwargs, :LBO, false)
    max_LBO_dim::Int = get(kwargs, :max_LBO_dim, 10)
    min_LBO_dim::Int = get(kwargs, :min_LBO_dim, 4)
  
    write_when_maxdim_exceeds::Union{Int,Nothing} = get(
      kwargs, :write_when_maxdim_exceeds, nothing
    )
  
    # eigsolve kwargs
    eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14)
    eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 3)
    eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 1)
    eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)
  
    # TODO: add support for non-Hermitian DMRG
    ishermitian::Bool = get(kwargs, :ishermitian, true)
  
    # TODO: add support for targeting other states with DMRG
    # (such as the state with the largest eigenvalue)
    # get(kwargs, :eigsolve_which_eigenvalue, :SR)
    eigsolve_which_eigenvalue::Symbol = :SR
  
    # TODO: use this as preferred syntax for passing arguments
    # to eigsolve
    #default_eigsolve_args = (tol = 1e-14, krylovdim = 3, maxiter = 1,
    #                         verbosity = 0, ishermitian = true,
    #                         which_eigenvalue = :SR)
    #eigsolve = get(kwargs, :eigsolve, default_eigsolve_args)
  
    # Keyword argument deprecations
    if haskey(kwargs, :maxiter)
      error("""maxiter keyword has been replaced by eigsolve_krylovdim.
               Note: compared to the C++ version of ITensor,
               setting eigsolve_krylovdim 3 is the same as setting
               a maxiter of 2.""")
    end
  
    if haskey(kwargs, :errgoal)
      error("errgoal keyword has been replaced by eigsolve_tol.")
    end
  
    if haskey(kwargs, :quiet)
      error("quiet keyword has been replaced by outputlevel")
    end
  
    psi = copy(psi0)
    N = length(psi)
  
    # Orthogonalize to site 1, which is the first site to be optimized
    if !isortho(psi) || ITensors.orthocenter(psi) != 1
      orthogonalize!(psi, 1)
    end
    @assert isortho(psi) && ITensors.orthocenter(psi) == 1
  
    # Position sets which indices of the MPO are hanging loose
    position!(PH, psi, 1)
    energy = 0.0
  
    # If doing LBO, initialize the set of Rs
    if LBO
      Rs = [ITensor() for _ in 1:N-1]
      PH_original = copy(PH)
      position!(PH_original, psi, 1)
    end
  
    # Iterate through the sweeps
    for sw in 1:nsweep(sweeps)
      sw_time = @elapsed begin
        maxtruncerr = 0.0
  
        if !isnothing(write_when_maxdim_exceeds) &&
          maxdim(sweeps, sw) > write_when_maxdim_exceeds
          if outputlevel >= 2
            println(
              "write_when_maxdim_exceeds = $write_when_maxdim_exceeds and maxdim(sweeps, sw) = $(maxdim(sweeps, sw)), writing environment tensors to disk",
            )
          end
          PH = disk(PH)
        end
  
        # Start iterating through the sites 
        #`b`is the bond number and `ha` is the half-sweep number
        # NOTE ha=1 on first half of the sweep, and ha=2 the other way
        for (b, ha) in sweepnext(N) 
          ITensors.@debug_check begin
            checkflux(psi)
            checkflux(PH)
          end
  
          """
          Optionally perform local basis optimization 
          """
  
          # NOTE: START WITH LARGER LBO_DIM and decrease with each sweep!
          if LBO && sw>1
  
            # Get the left indices of a tensor 
            function get_Linds(O::ITensor, tag::String)
              Rind = getfirst(x -> hastags(x, tag), inds(O))
              return noncommoninds([Rind], inds(O))
            end
  
            # Put into canonical form
            orthogonalize!(psi,b)
            M = psi[b] 
            H = PH_original.H[b] # we start fresh each time
  
            # Initialize the vector of Rs
            if !hastags(M,"Trunc")
              M,R = factorize(M, ortho="right", get_Linds(M, "Site"), tags="Trunc")
              Rs[b] = R
            else 
              R = Rs[b]
            end
  
            Linds = get_Linds(M, "Trunc")
  
            # Perform the first SVD
            X, Λ, Y = svd(M, Linds, lefttags="svd1,u", righttags="svd1,v")
            u = getfirst(x -> hastags(x, "svd1,u"), inds(Λ))
            # s = getfirst(x -> hastags(x, "Site"), inds(R))
            t = getfirst(x -> hastags(x, "Trunc"), inds(R))
  
            # Make and optimize R̃
            R̃ = R * Y * Λ 
            X̃, Λ̃, R = svd(R̃, u, lefttags="svd2,u", righttags="Trunc", maxdim=max_LBO_dim)
  
            t = getfirst(x -> hastags(x, "Trunc"), inds(R))
            # Make the new M̃ and H̃
            M̃ = Λ̃ * X̃ * X
            H̃ = dag(prime(R)) * H * R

            # Update
            psi[b] = M̃
            Rs[b] = R
            PH.H[b] = H̃
          end
          # Normalize? 
          normalize!(psi)
  
          @timeit_debug timer "dmrg: position!" begin
            # MPO hangs loose at site b now 
            position!(PH, psi, b)
          end
  
          ITensors.@debug_check begin
            checkflux(psi)
            checkflux(PH)
          end
  
          @timeit_debug timer "dmrg: psi[b]*psi[b+1]" begin
            # chunk together (contract) two sites 
            phi = psi[b] * psi[b + 1]
          end
  
          @timeit_debug timer "dmrg: eigsolve" begin
            # optimize the two-site chunk given the MPS
            vals, vecs = eigsolve(
              PH,
              phi,
              1,
              eigsolve_which_eigenvalue;
              ishermitian=ishermitian,
              tol=eigsolve_tol,
              krylovdim=eigsolve_krylovdim,
              maxiter=eigsolve_maxiter,
            )
          end
          energy::Number = vals[1]
          phi::ITensor = vecs[1]
  
          # If on first half of sweep, ortho centre is on the left
          # On second half of sweep, switch to right 
          ortho = ha == 1 ? "left" : "right"
  
          drho = nothing
          if noise(sweeps, sw) > 0.0
            @timeit_debug timer "dmrg: noiseterm" begin
              # Use noise term when determining new MPS basis
              drho = noise(sweeps, sw) * noiseterm(PH, phi, ortho)
            end
          end
  
          ITensors.@debug_check begin
            checkflux(phi)
          end
  
          # SVD the two-site chunk apart 
          # Optionally add a perturbation (noise)
          @timeit_debug timer "dmrg: replacebond!" begin
            spec = replacebond!(
              psi,
              b,
              phi;
              maxdim=maxdim(sweeps, sw),
              mindim=mindim(sweeps, sw),
              cutoff=cutoff(sweeps, sw),
              eigen_perturbation=drho,
              ortho=ortho,
              normalize=true,
              which_decomp=which_decomp,
              svd_alg=svd_alg,
            )
          end
          maxtruncerr = max(maxtruncerr, spec.truncerr)
  
          ITensors.@debug_check begin
            checkflux(psi)
            checkflux(PH)
          end
  
          if outputlevel >= 2
            @printf(
              "Sweep %d, half %d, bond (%d,%d) energy=%.12f\n", sw, ha, b, b + 1, energy
            )
            @printf(
              "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
              cutoff(sweeps, sw),
              maxdim(sweeps, sw),
              mindim(sweeps, sw)
            )
            @printf(
              "  Trunc. err=%.2E, bond dimension %d\n", spec.truncerr, dim(linkind(psi, b))
            )
            flush(stdout)
          end
  
          sweep_is_done = (b == 1 && ha == 2)
          measure!(
            obs;
            energy=energy,
            psi=psi,
            bond=b,
            sweep=sw,
            half_sweep=ha,
            spec=spec,
            outputlevel=outputlevel,
            sweep_is_done=sweep_is_done,
          )
        end
      end
      if outputlevel >= 1
        @printf(
          "After sweep %d energy=%.12f maxlinkdim=%d maxerr=%.2E time=%.3f\n",
          sw,
          energy,
          maxlinkdim(psi),
          maxtruncerr,
          sw_time
        )
        flush(stdout)
      end
      isdone = checkdone!(obs; energy=energy, psi=psi, sweep=sw, outputlevel=outputlevel)
  
      isdone && break
    end
  
    # Transform back to original basis if implementing LBO
    if LBO
        for b in 1:length(Rs)
          psi[b] *= Rs[b] # Transform back into original basis 
        end
        return energy, psi, Rs
    end
    return (energy, psi)
  end