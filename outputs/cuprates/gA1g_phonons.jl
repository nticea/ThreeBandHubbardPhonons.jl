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
Nx=16
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
ωA1=0.1
ωB1=0
gA1=0.05
gB1=0

## GLOBAL MODE CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
# Subtract 1 from this to get the maximum number of phonons allowed in that mode 
COPPER_DIM_1 = 1 
COPPER_DIM_2 = 1 
COPPER_DIM_3 = 1 
PX_DIM_1 = 3 
PX_DIM_2 = 1 
PX_DIM_3 = 1 
PY_DIM_1 = 1 
PY_DIM_2 = 1 
PY_DIM_3 = 1 

# DMRG parameters 
DMRG_numsweeps = 20 # total number of iterations 
DMRG_numsweeps_per_save = DMRG_numsweeps # Not saving, so it doesn't matter 
DMRG_maxdim = [200,400,600,800,1000,1200,1400,1600,1800,2000]
DMRG_cutoff = 1E-10
DMRG_LBO = false
max_lbo_dim = 12 

## SAVE OUT INFO ##
DMRG_numsweeps_per_save = 3 # If don't want to save regularly, just set this to DMRG_numsweeps
println("Running DMRG...")
dmrg_run(Nx, Ny, yperiodic, 
        μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, 
        ωB1, ωA1, gB1, gA1, 
        doping, 
        COPPER_DIM_1, 
        COPPER_DIM_2, 
        COPPER_DIM_3,
        PX_DIM_1, 
        PX_DIM_2, 
        PX_DIM_3,
        PY_DIM_1, 
        PY_DIM_2, 
        PY_DIM_3,
        DMRG_numsweeps, DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
        disk_save=false,
        dir_path=@__DIR__)

println("Computing correlations...")
correlations_run(Nx, Ny, yperiodic, μ, εd, εp, tpd, tpp, Vpd, Upp, Udd, 
                ωB1, ωA1, gB1, gA1, doping; dir_path=@__DIR__)