using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
include(joinpath(@__DIR__,"../src/model.jl"))
include(joinpath(@__DIR__,"../src/plotting.jl"))
include(joinpath(@__DIR__,"../src/utilities.jl"))
include(joinpath(@__DIR__,"../src/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/2-7-4_18:18:30_40Nx_1Ny.h5"

# If loading in minimal results, there are no states 
load_gs=false
load_phi=false
load_psi=false

println("Loading data...")
params, TBHModel, dmrg_results, eq_corr = load_structs(loadpath, 
                                                        load_gs=load_gs, 
                                                        load_phi=load_phi, 
                                                        load_psi=load_psi)

# Visualize things we want                                                        
plot_equilibrium_correlations(eq_corr, "spin")
plot_equilibrium_correlations(eq_corr, "charge")
plot_equilibrium_correlations(eq_corr, "sSC")
plot_equilibrium_correlations(eq_corr, "pSC")
plot_equilibrium_correlations(eq_corr, "dSC")
plot_charge_density(dmrg_results)
plot_phonon_density(dmrg_results)
plot_spin_density(dmrg_results)