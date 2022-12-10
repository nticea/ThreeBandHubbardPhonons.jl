using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/plotting.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))
include(joinpath(@__DIR__,"../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 

# gA1 
loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/cuprates/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0ωB1_0.2ωA1_0gB1_0.05gA1_results.h5"
# gB1 
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/cuprates/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0ωA1_0.005gB1_0gA1_results.h5"
# gA1 and gB1 
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/cuprates/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0.2ωA1_0.005gB1_0.05gA1_results.h5"

dmrg_results = load_dmrg_results_minimal(loadpath)
# Flags 
do_fit = true
do_save = true 

## PLOTTING ## 

# plot_charge_density(dmrg_results)
# plot_spin_density(dmrg_results)
# plot_phonon_density(dmrg_results) # plot sum of all nonzero modes 
# plot_densities(dmrg_results)

eq_corr = load_equilibrium_correlations(loadpath)
ploteq = plot_equilibrium_correlations(eq_corr, do_fit=do_fit)
plotd = plot_densities(dmrg_results)

if do_save 
    savefig(ploteq, "equilibrium_correlations.pdf")
    savefig(plotd, "densities.pdf")
end





