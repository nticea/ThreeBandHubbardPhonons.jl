function dmrg_run(Nx, Ny, yperiodic, μ, εd, εp, 
                tpd, tpp, Upd, Upp, Udd, ω, g0pp, g0dd, g1pd, 
                g1dp, g1pp, doping, max_phonons, DMRG_numsweeps,
                DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
                disk_save=false,
                dir_path=@__DIR__)

    param_stamp = "$(Nx)Nx_$(Ny)Ny_$(εp)εp_$(tpd)tpd_$(tpp)tpp_$(Upd)Upd_$(Upp)Upp_$(Udd)Udd_$(doping)doping_$(ω)ω_$(g0pp)g0pp_$(g0dd)g0dd_$(g1pd)g1pd_$(g1dp)g1dp_$(g1pp)g1pp"
    save_path = joinpath(dir_path,param_stamp*".h5")
    results_save_path = joinpath(@__DIR__,param_stamp*"_results.h5")

    # Create the output file 
    # output_path = joinpath(@__DIR__,param_stamp*"_out.log")
    # try
    #     open(output_path, "a") do s
    #         println(s, "Appending to file...")
    #     end
    # catch
    #     touch(output_path)
    #     open(output_path, "w") do s
    #         println(s, "Creating a new file...")
    #     end
    # end

    ## CODE ## 
    global params
    global dmrg_results 
    global TBHModel
    global eq_corr

    # Write everything to my output file 
    # open(output_path, "a") do io
    #     redirect_stdout(io) do 

    try
        global params = load_params(save_path)
        println("Loading parameters from ", save_path)
    catch 
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
    catch 
        println("Running DMRG...")
        global dmrg_results = run_DMRG(TBHModel, params, 
                                    DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, alg="divide_and_conquer", disk_save=disk_save)
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
                                    DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, alg="divide_and_conquer", disk_save=disk_save)
        dmrg_results_minimal = DMRGResultsMinimal(dmrg_results.ground_state_energy,
                                    dmrg_results.ground_state_entropy, 
                                    dmrg_results.charge_density,
                                    dmrg_results.phonon_density,
                                    dmrg_results.spin_density)
        save_structs(dmrg_results, save_path)
        save_structs(dmrg_results_minimal, results_save_path)
        println("Interim DMRG save")
    end
end
#     end
# end

function correlations_run(Nx, Ny, yperiodic, μ, εd, εp, 
                            tpd, tpp, Upd, Upp, Udd, ω, g0pp, g0dd, g1pd, 
                            g1dp, g1pp, doping, max_phonons, DMRG_numsweeps,
                            DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
                            dir_path=@__DIR__)

    param_stamp = "$(Nx)Nx_$(Ny)Ny_$(εp)εp_$(tpd)tpd_$(tpp)tpp_$(Upd)Upd_$(Upp)Upp_$(Udd)Udd_$(doping)doping_$(ω)ω_$(g0pp)g0pp_$(g0dd)g0dd_$(g1pd)g1pd_$(g1dp)g1dp_$(g1pp)g1pp"
    save_path = joinpath(dir_path,param_stamp*".h5")
    results_save_path = joinpath(@__DIR__,param_stamp*"_results.h5")

    # Create the output file 
    # output_path = joinpath(@__DIR__,param_stamp*"_out.log")
    # try
    #     open(output_path, "a") do s
    #         println(s, "Appending to file...")
    #     end
    # catch
    #     touch(output_path)
    #     open(output_path, "w") do s
    #         println(s, "Creating a new file...")
    #     end
    # end

    ## CODE ## 
    global eq_corr

    # Write everything to my output file 
    # open(output_path, "a") do io
    #     redirect_stdout(io) do 

    # Load in the parameters 
    params = load_params(save_path)
    println("Loading parameters from ", save_path)

    # Load the DMRG results         
    println("Loading DMRG results")
    dmrg_results = load_dmrg_results(save_path)

    # Load the mdoel 
    TBHModel = ThreeBandModel(params, dmrg_results) # load in the correct sites if we already have a wavefcn

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

    _,_,dSC_dpx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "d-px", "d-px", "dSC")
    global eq_corr.dSC_dpx = dSC_dpx 
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for d-px bond...")

    _,_,dSC_dydy = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dy-dy", "dy-dy", "dSC")
    global eq_corr.dSC_dydy = dSC_dydy 
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for dy-dy bond...")

    _,_,dSC_pyd = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-d", "py-d", "dSC")
    global eq_corr.dSC_pyd = dSC_pyd 
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py-d bond...")

    _,_,dSC_pypx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-px", "py-px", "dSC")
    global eq_corr.dSC_pypx = dSC_pypx 
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py-px bond...")

    _,_,dSC_py1px2 = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py1-px2", "py1-px2", "dSC")
    global eq_corr.dSC_py1px2 = dSC_py1px2 
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py1-px2 bond...")

    _,_,dSC_dxdx = compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dx-dx", "dx-dx", "dSC")
    global eq_corr.dSC_dxdx = dSC_dxdx 
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py-px bond...")

end