using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include(joinpath(@__DIR__, "../../src/model.jl"))
include(joinpath(@__DIR__, "../../src/plotting.jl"))
include(joinpath(@__DIR__, "../../src/utilities.jl"))
include(joinpath(@__DIR__, "../../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
lambda0_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/lambda_32Nx_2Ny/32Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0ωB1_0ωA1_0gB1_0gA1_results.h5"
lambda0_25_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/lambda_32Nx_2Ny/32Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1gB1_0gA1_results.h5"
lambda0_025_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/lambda_32Nx_2Ny/32Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.3162gB1_0gA1_results.h5"
lambda0_0025_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/lambda_32Nx_2Ny/32Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.1gB1_0gA1_results.h5"
lambda0_0025_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/lambda_32Nx_2Ny/32Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.03162gB1_0gA1_results.h5"

loadpath = lambda0_loadpath

# Flags 
do_fit = true
do_save = true

## PLOTTING ## 

dmrg_results = load_dmrg_results_minimal(loadpath)
plotd = plot_densities(dmrg_results)

try
    global eq_corr = load_equilibrium_correlations(loadpath)
    global ploteq = plot_equilibrium_correlations(eq_corr)
catch
    global ploteq = plot()
end

# λ = calculate_λ(params)
# @show λ

if do_save
    savefig(ploteq, "equilibrium_correlations.png")
    savefig(plotd, "densities.png")
end





