using ITensors.HDF5

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)

        # If the file already exists, then we either append to it or overwrite 
        if isfile(path)
            h5open(path, "r+") do file
                if haskey(file, n) #if key already exists, we want to rewrite 
                    delete_object(file, n)
                    write(file, n, d)
                else
                    write(file, n, d) 
                end
            end
        else # If the file does not exist, create it 
            h5open(path, "w") do file
                write(file, n, d) 
            end
        end
    end
end

function load_params(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return Parameters(d["N"], d["Nx"], d["Ny"], d["Nsites"], 
            d["yperiodic"], d["doping"],d["max_phonons"], d["init_phonons"],
            d["μ"], d["εd"],d["εp"], d["tpd"], d["tpp"], d["Vpd"],d["Upp"], d["Udd"],
            d["ω"], d["g0pp"],d["g0dd"], d["g1pd"], d["g1dp"], d["g1pp"],
            d["dim_copper_mode_1"], d["dim_copper_mode_2"], d["dim_copper_mode_3"],
            d["dim_oxygen_mode_1"], d["dim_oxygen_mode_2"], d["dim_oxygen_mode_3"],
            d["DMRG_numsweeps"], d["DMRG_noise"],d["DMRG_maxdim"],
            d["DMRG_cutoff"], d["DMRG_LBO"], d["max_LBO_dim"],d["min_LBO_dim"],
            d["mid"], d["T"], d["τ"], 
            d["TEBD_cutoff"], d["TEBD_maxdim"], d["TEBD_LBO"])
end

function load_dmrg_results(loadpath::String; load_gs::Bool=true)
    f = h5open(loadpath,"r")
    if load_gs
        ground_state = read(f, "ground_state", MPS)
    else
        ground_state = nothing
    end
    d = read(f)
    return DMRGResults(d["nsweep"], d["maxdim"], d["cutoff"],d["noise"], 
                            ground_state, d["ground_state_energy"], 
                            d["ground_state_entropy"], 0,
                            d["charge_density"], d["phonon_density"], d["spin_density"])
end

function load_equilibrium_correlations(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return EquilibriumCorrelations(d["start"], d["stop"], d["spin"], d["charge"], 
                                    d["dSC_dxdx"], d["dSC_dpx"], d["dSC_dydy"], 
                                    d["dSC_pyd"], d["dSC_pypx"], d["dSC_py1px2"])
end

function load_results(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    @show keys(d)
    dmrg_results = DMRGResultsMinimal(d["ground_state_energy"], 
                                    d["ground_state_entropy"], d["charge_density"], 
                                    d["phonon_density"], d["spin_density"])

    eq_corr = EquilibriumCorrelations(d["start"], d["stop"], d["spin"], d["charge"], 
                                    d["dSC_dxdx"], d["dSC_dpx"], d["dSC_dydy"], 
                                    d["dSC_pyd"], d["dSC_pypx"], d["dSC_py1px2"])
    return dmrg_results, eq_corr
end

function load_dmrg_results_minimal(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return DMRGResultsMinimal(d["ground_state_energy"], d["ground_state_entropy"], 
                            d["charge_density"], d["phonon_density"], d["spin_density"])
end