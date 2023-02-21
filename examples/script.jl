## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Dates
include(joinpath(@__DIR__, "../src/model.jl"))
include(joinpath(@__DIR__, "../src/utilities.jl"))
include(joinpath(@__DIR__, "../src/plotting.jl"))

ITensors.set_warn_order(50)

## PARAMETERS ## 

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

# Initialize 
println("Initializing...")
params = parameters(Nx=Nx, Ny=Ny, yperiodic=yperiodic, μ=μ, εd=εd, εp=εp, tpd=tpd,
    tpp=tpp, Vpd=Vpd, Upp=Upp, Udd=Udd,
    ωB1=ωB1, ωA1=ωA1, gB1=gB1, gA1=gA1, doping=doping,
    dim_copper_mode_1=COPPER_DIM_1, dim_copper_mode_2=COPPER_DIM_2,
    dim_copper_mode_3=COPPER_DIM_3,
    dim_oxygen_x_mode_1=PX_DIM_1, dim_oxygen_x_mode_2=PX_DIM_2,
    dim_oxygen_x_mode_3=PX_DIM_3, dim_oxygen_y_mode_1=PY_DIM_1,
    dim_oxygen_y_mode_2=PY_DIM_2, dim_oxygen_y_mode_3=PY_DIM_3,
    DMRG_numsweeps=DMRG_numsweeps, DMRG_noise=DMRG_noise,
    DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff);
# The Hamiltonian MPO 
TBHModel = ThreeBandModel(params)

# Run DMRG first few sweeps 
println("Finding ground state...")
global dmrg_results = run_DMRG(TBHModel, params, DMRG_numsweeps_per_save=DMRG_numsweeps_per_save,
    alg="divide_and_conquer", disk_save=false);

for i in 1:floor(Int, DMRG_numsweeps / DMRG_numsweeps_per_save)
    global dmrg_results = run_DMRG(dmrg_results, TBHModel, params,
        DMRG_numsweeps_per_save=DMRG_numsweeps_per_save,
        alg="divide_and_conquer", disk_save=false)
    @show sum(dmrg_results.phonon_density) / params.N
end
plot_densities(dmrg_results)

# Equilibrium correlations
println("Computing equilibrium correlations...")
start, stop, particle = @time compute_equilibrium_onsite_correlation(dmrg_results, TBHModel, params, "dx-dx", "particle")