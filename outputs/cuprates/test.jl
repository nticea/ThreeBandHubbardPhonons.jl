## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../../"))
include(joinpath(@__DIR__,"../../src/model.jl"))
include(joinpath(@__DIR__,"../../src/utilities.jl"))
include(joinpath(@__DIR__,"../../src/run.jl"))

# The equilibrium corrs require contracting tensors w many indices
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
ωB1g=0.1
ω1g=0.05
gB1g=0.01
gA1g=0.005

## GLOBAL MODE CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
# Subtract 1 from this to get the maximum number of phonons allowed in that mode 
COPPER_DIM_1 = 2 # maximum 1 phonon 
COPPER_DIM_2 = 2 # maximum 1 phonons 
COPPER_DIM_3 = 1 # maximum 0 phonons 
OXYGEN_DIM_1 = 2 # maximum 1 phonons
OXYGEN_DIM_2 = 2 # maximum 1 phonons 
OXYGEN_DIM_3 = 1 # maximum 1 phonons 

# DMRG parameters 
DMRG_numsweeps = 80 # total number of iterations 
DMRG_numsweeps_per_save = DMRG_numsweeps # Not saving, so it doesn't matter 
DMRG_maxdim = 2000
DMRG_cutoff = 1E-10
DMRG_LBO = false
max_lbo_dim = 12 

## SAVE OUT INFO ##
DMRG_numsweeps_per_save = 3 # If don't want to save regularly, just set this to DMRG_numsweeps
println("Running DMRG...")
dmrg_run(Nx, Ny, yperiodic, 
        μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, 
        ωB1g, ω1g, gB1g, gA1g, 
        doping, 
        COPPER_DIM_1, 
        COPPER_DIM_2, 
        COPPER_DIM_3,
        OXYGEN_DIM_1, 
        OXYGEN_DIM_2, 
        OXYGEN_DIM_3,
        DMRG_numsweeps, DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
        disk_save=true,
        dir_path=@__DIR__)

println("Computing correlations...")
correlations_run(Nx, Ny, yperiodic, μ, εd, εp, 
        tpd, tpp, Vpd, Upp, Udd, ω, g0pp, g0dd, g1pd, 
        g1dp, g1pp, doping, max_phonons, DMRG_numsweeps,
        DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
        dir_path=@__DIR__)