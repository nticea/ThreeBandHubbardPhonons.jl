## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using Dates
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))

## SAVING INFO ##
DO_SAVE = false
MINIMAL_SAVE = false

## PARAMETERS ## 

# Model 
Nx=4
Ny=2
yperiodic=true

μ=0
εd=0
εp=3
tpd=1
tpp=0.5
Upd=0
Upp=3
Udd=8
ω=0#0.5
g0pp=0#0.1
g0dd=0#0.1
g1pd=0#0.1
g1dp=0#0.1
g1pp=0#0.1
doping=0.125
max_phonons=0#1 # (n+1)*4 = total site dimension 

# DMRG parameters 
DMRG_numsweeps = 20
DMRG_maxdim = 2048
DMRG_cutoff = 1E-10 

## CODE ## 

# Save 
date_stamp = Dates.format(now(), "y-m-d_HH:MM:SS") 
param_stamp = "_$(Nx)Nx_$(Ny)Ny"
save_path_full = joinpath(@__DIR__,"../outputs",date_stamp*param_stamp*".h5")
save_path_minimal = joinpath(@__DIR__,"../outputs",date_stamp*param_stamp*"_minimal.h5")

if DO_SAVE
    println("Saving to ", save_path_full)
elseif MINIMAL_SAVE
    println("Saving to ", save_path_minimal)
end

# Initialize 
println("Initializing...")
params = parameters(Nx=Nx, Ny=Ny, yperiodic=yperiodic, μ=μ, εd=εd, εp=εp, tpd=tpd, tpp=tpp, Upd=Upd, 
                    Upp=Upp, Udd=Udd, ω=ω, g0pp=g0pp, g0dd=g0dd, g1pd=g1pd, 
                    g1dp=g1dp, g1pp=g1pp, doping=doping, 
                    max_phonons=max_phonons, DMRG_numsweeps=DMRG_numsweeps,
                    DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff)
TBHModel = ThreeBandModel(params)
if MINIMAL_SAVE
    save_structs(params, save_path_minimal)
end
if DO_SAVE
    save_structs(params, save_path_full)
    save_structs(TBHModel, save_path_full)
end

# Visualize the lattice
visualize_lattice(params)

# Run DMRG
println("Finding ground state...")
dmrg_results = run_DMRG(TBHModel, params, alg="divide_and_conquer")
if MINIMAL_SAVE
    save_structs(DMRGResultsMinimal(dmrg_results.ground_state_energy, 
                            dmrg_results.ground_state_entropy, 
                            dmrg_results.charge_density, dmrg_results.phonon_density,
                            dmrg_results.spin_density), save_path_minimal)
end
if DO_SAVE
    save_structs(dmrg_results, save_path_full)
end

# Equilibrium correlations
println("Computing equilibrium correlations...")
eq_corr = compute_all_equilibrium_correlations(dmrg_results, TBHModel, params)
if MINIMAL_SAVE
    save_structs(eq_corr, save_path_minimal)
end
if DO_SAVE
    save_structs(eq_corr, save_path_full)
end