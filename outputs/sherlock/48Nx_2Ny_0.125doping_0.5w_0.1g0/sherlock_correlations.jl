## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../../../"))
include(joinpath(@__DIR__,"../../../src/model.jl"))
include(joinpath(@__DIR__,"../../../src/utilities.jl"))

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
Upd=0
Upp=3
Udd=8
doping=0.125
ω=0.5
g0pp=0.1
g0dd=0.1
g1pd=0
g1dp=0
g1pp=0
max_phonons=4 # (n+1)*4 = total site dimension 

## SAVE OUT INFO ##
dir_path = "/scratch/users/nticea"
correlations_run(Nx=Nx, Ny=Ny, yperiodic=yperiodic, μ=μ, εd=εd, εp=εp, 
                tpd=tpd, tpp=tpp, Upd=Upd, 
                Upp=Upp, Udd=Udd, ω=ω, g0pp=g0pp, g0dd=g0dd, g1pd=g1pd, 
                g1dp=g1dp, g1pp=g1pp, doping=doping, 
                max_phonons=max_phonons, DMRG_numsweeps=DMRG_numsweeps,
                DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff,
                DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, dir_path=dir_path)
