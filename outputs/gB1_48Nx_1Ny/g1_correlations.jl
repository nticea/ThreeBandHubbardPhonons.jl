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
Nx = 48
Ny = 1
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
gB1 = 1


println("Computing correlations...")
correlations_run(Nx, Ny, yperiodic, μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
        ωB1, ωA1, gB1, gA1, doping; checkpoint_path="/scratch/users/nticea",
        results_path=@__DIR__)