using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/plotting.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))
include(joinpath(@__DIR__,"../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
loadpath_full = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/4Nx_2Ny_3εp_1tpd_0.5tpp_0Upd_3Upp_8Udd_0.125doping.h5"
loadpath_results = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/4Nx_2Ny_3εp_1tpd_0.5tpp_0Upd_3Upp_8Udd_0.125dopingresults.h5"
#loadpath = "/Users/nicole/sherlock/code/ThreeBandHubbardPhonons.jl/outputs/96Nx_2Ny_3εp_1tpd_0.5tpp_0Upd_3Upp_8Udd_0.125doping.h5"

println("Loading data...")
dmrg_results_minimal, eq_corr = load_results(loadpath_results)
dmrg_results = load_dmrg_results(loadpath_full)

# Visualize things we want                                                        
plot_equilibrium_correlations(eq_corr, "spin")
plot_equilibrium_correlations(eq_corr, "charge")
plot_equilibrium_correlations(eq_corr, "sSC")
plot_equilibrium_correlations(eq_corr, "pSC")
plot_equilibrium_correlations(eq_corr, "dSC")
plot_charge_density(dmrg_results)
plot_phonon_density(dmrg_results)
plot_spin_density(dmrg_results)