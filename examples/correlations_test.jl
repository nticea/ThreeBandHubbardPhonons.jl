## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Dates
include(joinpath(@__DIR__, "../src/model.jl"))
include(joinpath(@__DIR__, "../src/utilities.jl"))
include(joinpath(@__DIR__, "../src/plotting.jl"))
include(joinpath(@__DIR__, "../src/run.jl"))

# The equilibrium corrs require contracting tensors w many indices
ITensors.set_warn_order(50)

# Model 
Nx = 16
Ny = 1
yperiodic = false

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
gB1 = 0.1

λ = gB1^2 / (4 * ωB1)
@show λ

## GLOBAL MODE CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
# Subtract 1 from this to get the maximum number of phonons allowed in that mode 
COPPER_DIM_1 = 1
COPPER_DIM_2 = 1
COPPER_DIM_3 = 1
PX_DIM_1 = 1
PX_DIM_2 = 1
PX_DIM_3 = 1
PY_DIM_1 = 3
PY_DIM_2 = 1
PY_DIM_3 = 1

# DMRG parameters 
DMRG_numsweeps = 9 # total number of iterations 
DMRG_numsweeps_per_save = 3 # If don't want to save regularly, just set this to DMRG_numsweeps
DMRG_maxdim = [50, 50, 50, 50, 50,
    100, 100, 100, 100]
DMRG_noise = [1E-6, 1E-7, 1E-8, 1E-9]
DMRG_cutoff = 1E-12
overwrite_sweeps = false

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
    checkpoint_path=@__DIR__,
    results_path=@__DIR__)

println("Computing correlations...")
correlations_run_checkpoint(Nx, Ny, yperiodic, μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
    ωB1, ωA1, gB1, gA1, doping; checkpoint_path=@__DIR__,
    results_path=@__DIR__)

correlations_run(Nx, Ny, yperiodic, μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
    ωB1, ωA1, gB1, gA1, doping; checkpoint_path=@__DIR__,
    results_path=@__DIR__)

eq_corr_true = load_equilibrium_correlations("/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/examples/16Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.1gB1_0gA1_results.h5")
eq_corr_checkpoint = load_equilibrium_correlations("/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/examples/16Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.1gB1_0gA1_correlation_results.h5")

for fn in fieldnames(typeof(eq_corr_checkpoint))
    if getfield(eq_corr_true, fn) != getfield(eq_corr_checkpoint, fn)
        @show fn
    end
end