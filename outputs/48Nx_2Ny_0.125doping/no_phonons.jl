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
Nx=48
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
ω=0
g0pp=0
g0dd=0
g1pd=0
g1dp=0
g1pp=0
max_phonons=0 # (n+1)*4 = total site dimension 

# DMRG parameters 
DMRG_numsweeps = 80 # total number of iterations 
DMRG_maxdim = 2048
DMRG_cutoff = 1E-10

## SAVE OUT INFO ##
DMRG_numsweeps_per_save = 3 # If don't want to save regularly, just set this to DMRG_numsweeps
println("Running DMRG...")
dmrg_run(Nx, Ny, yperiodic, μ, εd, εp, 
        tpd, tpp, Vpd, Upp, Udd, ω, g0pp, g0dd, g1pd, 
        g1dp, g1pp, doping, max_phonons, DMRG_numsweeps,
        DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
        disk_save=true,
        dir_path=@__DIR__)

println("Computing correlations...")
correlations_run(Nx, Ny, yperiodic, μ, εd, εp, 
        tpd, tpp, Vpd, Upp, Udd, ω, g0pp, g0dd, g1pd, 
        g1dp, g1pp, doping, max_phonons, DMRG_numsweeps,
        DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
        dir_path=@__DIR__)