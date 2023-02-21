using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
include(joinpath(@__DIR__, "../../src/model.jl"))
include(joinpath(@__DIR__, "../../src/plotting.jl"))
include(joinpath(@__DIR__, "../../src/utilities.jl"))
include(joinpath(@__DIR__, "../../src/sites/site_hubbholst.jl"))
using ITensors
using ITensors.HDF5

# Load the results in 
g0_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_16Nx_1Ny/16Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_0ωB1_0ωA1_0gB1_0gA1_results.h5"
g1_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_16Nx_1Ny/16Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1gB1_0gA1_results.h5"
g15_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_16Nx_1Ny/16Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1.5gB1_0gA1_results.h5"
g175_loadpath = "/Users/nicole/Dropbox/Grad school/Devereaux lab/Hubbard model/3BHPhonons/ThreeBandHubbardPhonons.jl/outputs/gB1_16Nx_1Ny/16Nx_1Ny_3εp_1tpd_0.5tpp_0Vpd_3Upp_8Udd_0.125doping_1ωB1_0ωA1_1.75gB1_0gA1_results.h5"

# Flags 
do_fit = true
do_save = true

## PLOTTING ## 

dmrg_results = load_dmrg_results_minimal(loadpath)
eq_corr = load_equilibrium_correlations(loadpath)

cmap = cgrad(:Set1_4, 4, categorical=true)
gs = [0, 1, 1.5, 1.75]
loadpaths = [g0_loadpath, g1_loadpath, g15_loadpath, g175_loadpath]

global p = plot()
for (i, (g, loadpath)) in enumerate(zip(gs, loadpaths))
    corrs = load_equilibrium_correlations(loadpath).spin
    corrs = corrs[1:end] # discard the correlation at zero distance (for now...)
    xrange = collect(1:length(corrs))

    # Do fits
    type, a, b, fit = compare_fits(xrange, abs.(corrs))
    global p = plot!(p, log10.(xrange), log10.(abs.(corrs)), label="gB1=$(g)", c=cmap[i])
    global p = scatter(p, log10.(xrange), log10.(abs.(corrs)), label=nothing)
    global p = plot!(p, log10.(xrange), log10.(fit), label="$type, k=$b", c=cmap[i])
end

global p = title!(p, "16x1 chain, spin correlation")
plot(p)





