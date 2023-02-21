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
ωB1 = 0
gA1 = 0
gB1 = 0

## GLOBAL MODE CONSTANTS -- CAN'T THINK OF A GOOD WAY TO INCORPORATE THEM OTHERWISE!! ## 
# Subtract 1 from this to get the maximum number of phonons allowed in that mode 
COPPER_DIM_1 = 1
COPPER_DIM_2 = 1
COPPER_DIM_3 = 1
PX_DIM_1 = 1
PX_DIM_2 = 1
PX_DIM_3 = 1
PY_DIM_1 = 1
PY_DIM_2 = 1
PY_DIM_3 = 1

println("Computing correlations...")
correlations_run(Nx, Ny, yperiodic, μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
        ωB1, ωA1, gB1, gA1, doping; checkpoint_path="/scratch/users/nticea",
        results_path=@__DIR__)