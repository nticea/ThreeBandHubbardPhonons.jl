using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/plotting.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))
include(joinpath(@__DIR__,"../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
#loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/cuprates/48Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0.05ωA1_0gB1_0.005gA1_copy_results.h5"
loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/cuprates/48Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0.05ωA1_0.01gB1_0.005gA1_copy_results.h5"

println("Loading data...")
dmrg_results = load_dmrg_results_minimal(loadpath)

plot_charge_density(dmrg_results)
plot_spin_density(dmrg_results)
plot_phonon_density(dmrg_results, 1) # plot mode 1 
# plot_phonon_density(dmrg_results, 2) # plot mode 2 
# plot_phonon_density(dmrg_results, 3) # plot mode 3 

# do_fit = true
# eq_corr = load_equilibrium_correlations(loadpath)
# plot_equilibrium_correlations(eq_corr, "spin", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "charge", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "dSC_dxdx", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "dSC_dpx", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "dSC_dydy", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "dSC_pyd", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "dSC_pypx", do_fit=do_fit)
# plot_equilibrium_correlations(eq_corr, "dSC_py1px2", do_fit=do_fit)






