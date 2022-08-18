## IMPORTS ##
using Pkg
Pkg.activate(joinpath(@__DIR__,"../../../"))
using Dates
include(joinpath(@__DIR__,"../../../src/model.jl"))
include(joinpath(@__DIR__,"../../../src/utilities.jl"))

## PARAMETERS ## 

# Model 
Nx=96
Ny=2
yperiodic=true

μ=0
εd=0
εp=3
tpd=-1
tpp=-0.5
Upd=0
Upp=3
Udd=8
doping=0.125
ω=0
g0pp=0
g0dd=0
g1pd=0
g1dp=0
g1pp=0
max_phonons=0 # (n+1)*4 = total site dimension 

# DMRG parameters 
DMRG_numsweeps = 60 # total number of iterations 
DMRG_maxdim = 2048
DMRG_cutoff = 1E-10

## SAVE OUT INFO ##
DMRG_numsweeps_per_save = 5 # If don't want to save regularly, just set this to DMRG_numsweeps
param_stamp = "$(Nx)Nx_$(Ny)Ny_$(εp)εp_$(tpd)tpd_$(tpp)tpp_$(Upd)Upd_$(Upp)Upp_$(Udd)Udd_$(doping)doping_$(ω)ω_$(g0pp)g0pp_$(g0dd)g0dd_$(g1pd)g1pd_$(g1dp)g1dp_$(g1pp)g1pp"
save_path = joinpath("/scratch/users/nticea",param_stamp*".h5") #joinpath(@__DIR__,param_stamp*".h5")
results_save_path = joinpath(@__DIR__,param_stamp*"_results.h5")
output_path = joinpath(@__DIR__,param_stamp*"_out.log")

# Create the output file 
try
    open(output_path, "a") do s
        println(s, "Appending to file...")
    end
catch
    touch(output_path)
    open(output_path, "w") do s
        println(s, "Creating a new file...")
    end
end

## CODE ## 
global params
global dmrg_results 
global TBHModel
global eq_corr

# Write everything to my output file 
open(output_path, "a") do io
    redirect_stdout(io) do 

        try
            global params = load_params(save_path)
            println("Loading parameters from ", save_path)
        catch e
            @show e
            global params = parameters(Nx=Nx, Ny=Ny, yperiodic=yperiodic, μ=μ, εd=εd, εp=εp, 
                tpd=tpd, tpp=tpp, Upd=Upd, 
                Upp=Upp, Udd=Udd, ω=ω, g0pp=g0pp, g0dd=g0dd, g1pd=g1pd, 
                g1dp=g1dp, g1pp=g1pp, doping=doping, 
                max_phonons=max_phonons, DMRG_numsweeps=DMRG_numsweeps,
                DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff)
            println("Initializing parameters...")
            save_structs(params, save_path)
            println("Saving parameters to ", save_path)
        end

        # Initialize the model (sites and MPO)
        global TBHModel = ThreeBandModel(params)

        # Do the first DMRG run        
        try 
            println("Loading DMRG results")
            global dmrg_results = load_dmrg_results(save_path)
            global TBHModel = ThreeBandModel(params, dmrg_results) # load in the correct sites if we already have a wavefcn
        catch e
            @show e
            println("Running DMRG...")
            global dmrg_results = run_DMRG(TBHModel, params, 
                                        DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, alg="divide_and_conquer")
            dmrg_results_minimal = DMRGResultsMinimal(dmrg_results.ground_state_energy,
                                        dmrg_results.ground_state_entropy, 
                                        dmrg_results.charge_density,
                                        dmrg_results.phonon_density,
                                        dmrg_results.spin_density)
            save_structs(dmrg_results, save_path)
            save_structs(dmrg_results_minimal, results_save_path)
            println("Interim DMRG save")
        end

        # Do the rest of the DMRG runs 
        for _ in 1:floor(Int, DMRG_numsweeps/DMRG_numsweeps_per_save)
            global dmrg_results = run_DMRG(dmrg_results, TBHModel, params, 
                                        DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, alg="divide_and_conquer")
            dmrg_results_minimal = DMRGResultsMinimal(dmrg_results.ground_state_energy,
                                        dmrg_results.ground_state_entropy, 
                                        dmrg_results.charge_density,
                                        dmrg_results.phonon_density,
                                        dmrg_results.spin_density)
            save_structs(dmrg_results, save_path)
            save_structs(dmrg_results_minimal, results_save_path)
            println("Interim DMRG save")
        end

        # Compute and save the correlations 
        try 
            global eq_corr = load_equilibrium_correlations(results_save_path)
            println("Loading equilibrium correlations")
        catch e 
            @show e   
            global eq_corr = EquilibriumCorrelations(0,0,[0],[0],[0],[0],[0])
            println("Computing equilibrium correlations")
        end

        spin_corr, start, stop  = compute_equilibrium_correlation(dmrg_results, TBHModel, params, corrtype="spin")
        global eq_corr.start = start
        global eq_corr.stop = stop 
        global eq_corr.spin = spin_corr 
        save_structs(eq_corr, results_save_path)
        println("Saving spin correlation...")

        charge_corr,_,_ = compute_equilibrium_correlation(dmrg_results, TBHModel, params, corrtype="charge")
        global eq_corr.charge = charge_corr 
        save_structs(eq_corr, results_save_path)
        println("Saving charge correlation...")

        sSC_corr,_,_ = compute_equilibrium_correlation(dmrg_results, TBHModel, params, corrtype="sSC")
        global eq_corr.sSC = sSC_corr 
        save_structs(eq_corr, results_save_path)
        println("Saving sSC correlation...")

        pSC_corr,_,_ = compute_equilibrium_correlation(dmrg_results, TBHModel, params, corrtype="pSC")
        global eq_corr.pSC = pSC_corr 
        save_structs(eq_corr, results_save_path)
        println("Saving pSC correlation...")

        dSC_corr,_,_ = compute_equilibrium_correlation(dmrg_results, TBHModel, params, corrtype="dSC")
        global eq_corr.dSC = dSC_corr 
        save_structs(eq_corr, results_save_path)
        println("Saving dSC correlation...")

    end
end