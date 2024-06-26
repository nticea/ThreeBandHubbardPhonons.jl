using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include(joinpath(@__DIR__, "../../src/model.jl"))
include(joinpath(@__DIR__, "../../src/plotting.jl"))
include(joinpath(@__DIR__, "../../src/utilities.jl"))
include(joinpath(@__DIR__, "../../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
g0_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_1Ny/48Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0ωB1_0ωA1_0gB1_0gA1_results.h5"
g1_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_1Ny/48Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1gB1_0gA1_results.h5"
g149_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_1Ny/48Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1.49gB1_0gA1_results.h5"
g15_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_1Ny/48Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1.5gB1_0gA1_results.h5"

# Flags 
do_fit = true
do_save = true

## PLOTTING ## 

gs = [0, 1, 1.49]
loadpaths = [g0_loadpath, g1_loadpath, g149_loadpath]

plot_multiple_densities(loadpaths, gs)





