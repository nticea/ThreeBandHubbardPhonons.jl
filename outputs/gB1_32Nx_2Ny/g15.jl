## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
include(joinpath(@__DIR__, "../../src/model.jl"))
include(joinpath(@__DIR__, "../../src/utilities.jl"))
include(joinpath(@__DIR__, "../../src/run.jl"))

# The equilibrium corrs require contracting tensors w many indices
ITensors.set_warn_order(50)

## PARAMETERS ## 

# Model 
Nx = 32
Ny = 2
yperiodic = true

μ = 0
εd = 0
εp = 3
tpd = 1
tpp = 0.5
Vpd = 0
Upp = 3
Udd = 8
doping = 0.125
ωA1 = 0
ωB1 = 1
gA1 = 0
gB1 = 1.5

## GLOBAL MODE CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
# Subtract 1 from this to get the maximum number of phonons allowed in that mode 
COPPER_DIM_1 = 1
COPPER_DIM_2 = 1
COPPER_DIM_3 = 1
PX_DIM_1 = 1
PX_DIM_2 = 1
PX_DIM_3 = 1
PY_DIM_1 = 15
PY_DIM_2 = 1
PY_DIM_3 = 1

# DMRG parameters 
DMRG_numsweeps = 80 # total number of iterations 
DMRG_numsweeps_per_save = 1 # If don't want to save regularly, just set this to DMRG_numsweeps
DMRG_maxdim = [50, 50, 50, 50, 50,
        100, 100, 100, 100, 100,
        200, 200, 200, 200, 200,
        300, 300, 300, 300, 300,
        500, 500, 500, 500, 500,
        700, 700, 700, 700, 700,
        900, 900, 900, 900, 900,
        1000, 1000, 1000, 1000, 1000,
        1500, 1500, 1500, 1500, 1500,
        1500, 1500, 1500, 1500, 1500,
        2000, 2000, 2000, 2000, 2000,
        2000, 2000, 2000, 2000, 2000,
        2500, 2500, 2500, 2500, 2500,
        2500, 2500, 2500, 2500, 2500,
        3000, 3000, 3000, 3000, 3000]
DMRG_noise = [1E-6, 1E-7, 1E-8, 1E-9, 0,
        1E-6, 1E-7, 1E-8, 1E-9, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        0, 0, 0, 0, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        0, 0, 0, 0, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        0, 0, 0, 0, 0,
        1E-7, 1E-8, 1E-9, 1E-10, 0,
        0, 0, 0, 0, 0]
DMRG_cutoff = 1E-12
overwrite_sweeps = false

# DMRG_numsweeps = 30 # total number of iterations 
# DMRG_numsweeps_per_save = 3 # If don't want to save regularly, just set this to DMRG_numsweeps
# DMRG_maxdim = [4000]
# DMRG_noise = [1E-7, 1E-8, 1E-9, 1E-10, 0]
# DMRG_cutoff = 1E-12
# overwrite_sweeps = true

# DMRG_numsweeps = 20 # total number of iterations 
# DMRG_numsweeps_per_save = 3
# DMRG_maxdim = [2500]
# DMRG_noise = [1E-6, 1E-7, 1E-8, 1E-9, 0]
# DMRG_cutoff = 1E-12

## SAVE OUT INFO ##
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
        DMRG_numsweeps, DMRG_noise,
        DMRG_maxdim, DMRG_cutoff,
        DMRG_numsweeps_per_save;
        overwrite_sweeps=overwrite_sweeps,
        disk_save=false,
        checkpoint_path="/scratch/users/nticea",
        results_path=@__DIR__)

println("Computing correlations...")
correlations_run(Nx, Ny, yperiodic, μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
        ωB1, ωA1, gB1, gA1, doping; checkpoint_path="/scratch/users/nticea",
        results_path=@__DIR__)