using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
include(joinpath(@__DIR__,"../../src/model.jl"))
include(joinpath(@__DIR__,"../../src/plotting.jl"))
include(joinpath(@__DIR__,"../../src/utilities.jl"))
include(joinpath(@__DIR__,"../../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 

# no phonons 
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_2Ny/48Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0ωB1_0ωA1_0gB1_0gA1_results.h5"
# 0.25
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_2Ny/48Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1gB1_0gA1_results.h5"
# 0.025 
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_2Ny/48Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0ωA1_0.1gB1_0gA1_results.h5"
# 0.0025
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.1gB1_0gA1_results.h5"
# 0.00025
loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_48Nx_2Ny/48Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0ωA1_0.01gB1_0gA1_results.h5"


# Flags 
do_fit = true
do_save = true 
do_eq_corr = false

## PLOTTING ## 

plotd = plot_densities(dmrg_results)
dmrg_results = load_dmrg_results_minimal(loadpath)

if do_eq_corr
    eq_corr = load_equilibrium_correlations(loadpath)
    ploteq = plot_equilibrium_correlations(eq_corr)
end

# λ = calculate_λ(params)
# @show λ

if do_save
    if do_eq_corr 
        savefig(ploteq, "equilibrium_correlations.png")
    end
    savefig(plotd, "densities.png")
end





