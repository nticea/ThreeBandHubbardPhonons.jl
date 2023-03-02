function dmrg_run(Nx, Ny, yperiodic,
    μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
    ωB1, ωA1, gB1, gA1,
    doping,
    dim_copper_mode_1,
    dim_copper_mode_2,
    dim_copper_mode_3,
    dim_oxygen_x_mode_1,
    dim_oxygen_x_mode_2,
    dim_oxygen_x_mode_3,
    dim_oxygen_y_mode_1,
    dim_oxygen_y_mode_2,
    dim_oxygen_y_mode_3,
    DMRG_numsweeps, DMRG_noise,
    DMRG_maxdim, DMRG_cutoff, DMRG_numsweeps_per_save;
    overwrite_sweeps=false,
    disk_save=false,
    checkpoint_path=@__DIR__,
    results_path=@__DIR__)

    ## CODE ## 
    global dmrg_results
    global TBHModel
    global eq_corr

    # Re-initializing the parameters every time in case we've updated the DMRG parameters
    params = parameters(Nx=Nx, Ny=Ny, yperiodic=yperiodic, μ=μ, εd=εd, εp=εp, tpd=tpd,
        tpp=tpp, Vpd=Vpd, Upp=Upp, Udd=Udd,
        ωB1=ωB1, ωA1=ωA1, gB1=gB1, gA1=gA1, doping=doping,
        dim_copper_mode_1=dim_copper_mode_1, dim_copper_mode_2=dim_copper_mode_2,
        dim_copper_mode_3=dim_copper_mode_3,
        dim_oxygen_x_mode_1=dim_oxygen_x_mode_1, dim_oxygen_x_mode_2=dim_oxygen_x_mode_2,
        dim_oxygen_x_mode_3=dim_oxygen_x_mode_3,
        dim_oxygen_y_mode_1=dim_oxygen_y_mode_1, dim_oxygen_y_mode_2=dim_oxygen_y_mode_2,
        dim_oxygen_y_mode_3=dim_oxygen_y_mode_3,
        DMRG_numsweeps=DMRG_numsweeps, DMRG_noise=DMRG_noise,
        DMRG_maxdim=DMRG_maxdim, DMRG_cutoff=DMRG_cutoff)
    println("Initializing parameters...")

    λ = calculate_λ(params)
    @show λ

    ## SAVING ## 
    param_stamp = "$(Nx)Nx_$(Ny)Ny_$(εp)εp_$(tpd)tpd_$(tpp)tpp_$(Vpd)Vpd_$(Upp)Upp_$(Udd)Udd_$(doping)doping_$(ωB1)ωB1_$(ωA1)ωA1_$(gB1)gB1_$(gA1)gA1"
    save_path = joinpath(checkpoint_path, param_stamp * ".h5")
    results_save_path = joinpath(results_path, param_stamp * "_results.h5")

    save_structs(params, save_path)
    save_structs(params, results_save_path)

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
    for _ in 1:floor(Int, DMRG_numsweeps / DMRG_numsweeps_per_save)
        global dmrg_results = run_DMRG(dmrg_results, TBHModel, params,
            DMRG_numsweeps_per_save=DMRG_numsweeps_per_save, alg="divide_and_conquer",
            overwrite_sweeps=overwrite_sweeps, disk_save=disk_save)
        # benchmarking 
        ψ_gs = dmrg_results.ground_state
        mldim = maxlinkdim(ψ_gs)

        dmrg_results_minimal = DMRGResultsMinimal(dmrg_results.ground_state_energy,
            dmrg_results.ground_state_entropy,
            dmrg_results.charge_density,
            dmrg_results.phonon_density,
            dmrg_results.spin_density)
        save_structs(dmrg_results, save_path)
        save_structs(dmrg_results_minimal, results_save_path)

        # save also an extra file for each maxdim
        maxdim_results_save_path = joinpath(results_path, param_stamp * "_$(mldim)maxdim" * "_results.h5")
        save_structs(dmrg_results_minimal, maxdim_results_save_path)

        println("Interim DMRG save")
    end
end

function correlations_run(Nx, Ny, yperiodic,
    μ, εd, εp, tpd, tpp, Vpd, Upp, Udd,
    ωB1, ωA1, gB1, gA1,
    doping;
    checkpoint_path=@__DIR__,
    results_path=@__DIR__)

    param_stamp = "$(Nx)Nx_$(Ny)Ny_$(εp)εp_$(tpd)tpd_$(tpp)tpp_$(Vpd)Vpd_$(Upp)Upp_$(Udd)Udd_$(doping)doping_$(ωB1)ωB1_$(ωA1)ωA1_$(gB1)gB1_$(gA1)gA1"
    save_path = joinpath(checkpoint_path, param_stamp * ".h5")
    results_save_path = joinpath(results_path, param_stamp * "_results.h5")

    ## CODE ## 
    global eq_corr

    # Load in the parameters 
    p = load_params(save_path)
    println("Loading parameters from ", save_path)

    # Load the DMRG results         
    println("Loading DMRG results")
    dmrg_results = load_dmrg_results(save_path)

    # Load the mdoel 
    HM = ThreeBandModel(p, dmrg_results) # load in the correct sites if we already have a wavefcn

    # Compute and save the correlations 
    try
        global eq_corr = load_equilibrium_correlations(results_save_path)
        println("Loading equilibrium correlations")
    catch e
        @show e
        global eq_corr = EquilibriumCorrelations(0, 0, [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0])
        println("Computing equilibrium correlations")
    end
    flush(stdout)

    start, stop, dSC_dxdx = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dx-dx", "dx-dx", "dSC")
    global eq_corr.start = start
    global eq_corr.stop = stop
    global eq_corr.dSC_dxdx = dSC_dxdx
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for dx-dx bond...")
    flush(stdout)

    _, _, spin_corr = @time compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "spin")
    global eq_corr.spin = spin_corr
    save_structs(eq_corr, results_save_path)
    println("Saving spin correlation...")
    flush(stdout)

    _, _, charge_corr = @time compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "charge")
    global eq_corr.charge = charge_corr
    save_structs(eq_corr, results_save_path)
    println("Saving charge correlation...")
    flush(stdout)

    _, _, particle_corr = @time compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "particle")
    global eq_corr.particle = particle_corr
    save_structs(eq_corr, results_save_path)
    println("Saving particle-particle correlation...")
    flush(stdout)

    _, _, sSC_corr = @time compute_equilibrium_onsite_correlation(dmrg_results, HM, p, "dx-dx", "sSC")
    global eq_corr.sSC = sSC_corr
    save_structs(eq_corr, results_save_path)
    println("Saving sSC-sSC correlation...")
    flush(stdout)

    _, _, pSC_dxdx = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dx-dx", "dx-dx", "pSC")
    global eq_corr.pSC_dxdx = pSC_dxdx
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for dy-dy bond...")
    flush(stdout)

    _, _, dSC_dpx = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "d-px", "d-px", "dSC")
    global eq_corr.dSC_dpx = dSC_dpx
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for d-px bond...")
    flush(stdout)

    _, _, dSC_pyd = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-d", "py-d", "dSC")
    global eq_corr.dSC_pyd = dSC_pyd
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py-d bond...")
    flush(stdout)

    _, _, dSC_pypx = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py-px", "py-px", "dSC")
    global eq_corr.dSC_pypx = dSC_pypx
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py-px bond...")
    flush(stdout)

    _, _, dSC_py1px2 = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "py1-px2", "py1-px2", "dSC")
    global eq_corr.dSC_py1px2 = dSC_py1px2
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for py1-px2 bond...")
    flush(stdout)

    _, _, dSC_dydy = @time compute_equilibrium_pairfield_correlation(dmrg_results, HM, p, "dy-dy", "dy-dy", "dSC")
    global eq_corr.dSC_dydy = dSC_dydy
    save_structs(eq_corr, results_save_path)
    println("Saving dSC correlations for dy-dy bond...")
    flush(stdout)
end