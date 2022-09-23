using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/plotting.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))
include(joinpath(@__DIR__,"../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
#loadpath = "/Users/nicole/sherlock/code/ThreeBandHubbardPhonons.jl/outputs/sherlock/48Nx_2Ny_0.125doping_0.5w_0.1g0/48Nx_2Ny_3εp_1tpd_0.5tpp_0Upd_3Upp_8Udd_0.125doping_0.5ω_0.1g0pp_0.1g0dd_0g1pd_0g1dp_0g1pp_results.h5"
loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/sherlock/32Nx_2Ny_0.125doping/8Nx_2Ny_3εp_1tpd_0.5tpp_0Upd_3Upp_8Udd_0.125doping_0ω_0g0pp_0g0dd_0g1pd_0g1dp_0g1pp_results.h5"

println("Loading data...")
dmrg_results = load_dmrg_results_minimal(loadpath)
eq_corr = load_equilibrium_correlations(loadpath)

# Visualize things we want
println("Plotting...")                                                        
plot_equilibrium_correlations(eq_corr, "spin")
plot_equilibrium_correlations(eq_corr, "charge")
plot_equilibrium_correlations(eq_corr, "dSC_dxdx")
plot_equilibrium_correlations(eq_corr, "dSC_dpx")
plot_equilibrium_correlations(eq_corr, "dSC_dydy")
plot_equilibrium_correlations(eq_corr, "dSC_pyd")
plot_equilibrium_correlations(eq_corr, "dSC_pypx")
plot_equilibrium_correlations(eq_corr, "dSC_py1px2")

plot_charge_density(dmrg_results)
plot_phonon_density(dmrg_results)
plot_spin_density(dmrg_results)

