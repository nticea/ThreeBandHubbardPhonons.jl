using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include(joinpath(@__DIR__, "../../src/model.jl"))
include(joinpath(@__DIR__, "../../src/plotting.jl"))
include(joinpath(@__DIR__, "../../src/utilities.jl"))
include(joinpath(@__DIR__, "../../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 

# no phonons 
loadpath_0 = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0ωB1_0ωA1_0gB1_0gA1_results.h5"
# 0.25
loadpath_0_25 = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1gB1_0gA1_results.h5"
# 0.025 
loadpath_0_025 = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0ωA1_0.1gB1_0gA1_results.h5"
# 0.0025
loadpath_0_0025 = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_0.1gB1_0gA1_results.h5"
# 0.00025
loadpath_0_00025 = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1/16Nx_2Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0.1ωB1_0ωA1_0.01gB1_0gA1_results.h5"

# Flags 
do_fit = true
do_save = true
do_eq_corr = true

## LOAD RESULTS AND PLOT ##
dmrg_results = load_dmrg_results_minimal(loadpath_0_00025)
plotd = plot_densities(dmrg_results)

cmap = cgrad(:Set1_5, 5, categorical=true)
λs = [0, 0.25, 0.025, 0.025, 0.00025]
loadpaths = [loadpath_0, loadpath_0_25, loadpath_0_025, loadpath_0_0025, loadpath_0_00025]

global p = plot()
for (i, (λ, loadpath)) in enumerate(zip(λs, loadpaths))
    corrs = load_equilibrium_correlations(loadpath).dSC_py1px2
    corrs = corrs[1:end] # discard the correlation at zero distance (for now...)
    xrange = collect(1:length(corrs))

    # Do fits
    type, a, b, fit = compare_fits(xrange, abs.(corrs))
    global p = plot!(p, log10.(xrange), log10.(abs.(corrs)), label="λ=$(λ)", c=cmap[i])
    global p = scatter(p, log10.(xrange), log10.(abs.(corrs)), label=nothing)
    global p = plot!(p, log10.(xrange), log10.(fit), label="$type, k=$b", c=cmap[i])
end

global p = title!(p, "16x2 chain, dSC_py1px2 correlation")
plot(p)

savefig(p, "16x2_dSC_py1px2.png")






# ## SAVING ## 
# if do_save
#     if do_eq_corr
#         savefig(ploteq, "equilibrium_correlations.png")
#     end
#     savefig(plotd, "densities.png")
# end





