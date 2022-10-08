## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using Dates
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))

ITensors.set_warn_order(50)

## PARAMETERS ## 

# Model 
Nx=8
Ny=2
yperiodic=true

μ=0
εd=0
εp=3
tpd=1
tpp=0.5
Vpd=0
Upp=3
Udd=8
doping=0.125
ω=0.1
g0pp=0.01
g0dd=0.01
g1pd=0.001
g1dp=0.001
g1pp=0.001
max_phonons=3 # (n+1)*4 = total site dimension 

## GLOBAL MODE CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
COPPER_DIM_1 = 1 # maximum 1 phonon 
COPPER_DIM_2 = 1 # maximum 2 phonon 
COPPER_DIM_3 = 1 # maximum 1 phonon 
OXYGEN_DIM_1 = 1 # maximum 1 phonon 
OXYGEN_DIM_2 = 1 # maximum 1 phonon 
OXYGEN_DIM_3 = 1 # maximum 1 phonon 

# DMRG parameters 
DMRG_numsweeps = 15 # total number of iterations 
DMRG_numsweeps_per_save = DMRG_numsweeps # Not saving, so it doesn't matter 
DMRG_maxdim = 64
DMRG_cutoff = 1E-10
DMRG_LBO = false
max_lbo_dim = 12 

# Initialize 
println("Initializing...")
params = parameters(Nx=Nx, Ny=Ny, yperiodic=yperiodic, μ=μ, εd=εd, εp=εp, tpd=tpd, tpp=tpp, Vpd=Vpd, 
                    Upp=Upp, Udd=Udd, ω=ω, g0pp=g0pp, g0dd=g0dd, g1pd=g1pd, 
                    g1dp=g1dp, g1pp=g1pp, doping=doping, 
                    dim_copper_mode_1=COPPER_DIM_1, dim_copper_mode_2=COPPER_DIM_2, 
                    dim_copper_mode_3=COPPER_DIM_3,
                    dim_oxygen_mode_1=OXYGEN_DIM_1, dim_oxygen_mode_2=OXYGEN_DIM_2, 
                    dim_oxygen_mode_3=OXYGEN_DIM_3,
                    DMRG_numsweeps=DMRG_numsweeps,
                    DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff, DMRG_LBO=DMRG_LBO);
# The Hamiltonian MPO 
TBHModel = ThreeBandModel(params);

# Run DMRG
println("Finding ground state...")
dmrg_results = run_DMRG(TBHModel, params, DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, 
                        alg="divide_and_conquer", disk_save=true);

# Equilibrium correlations
println("Computing equilibrium correlations...")
eq_corr = compute_all_equilibrium_correlations(dmrg_results, TBHModel, params)