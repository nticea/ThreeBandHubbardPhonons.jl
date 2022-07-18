using ITensors.HDF5

function save_structs(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)

        ## DEBUGGING ## 
        if n=="optimized_basis" && !isnothing(d)
            try 
                h5open(path, "r+") do file
                    for i in 1:length(d)
                        write(file, n*"_"*string(i), d[i])
                    end
                end
            catch e
                if !isfile(path)
                    h5open(path, "w") do file
                        for i in 1:length(d)
                            write(file, n*"_"*string(i), d[i])
                        end
                    end
                else
                    h5open(path, "r+") do file
                        for i in 1:length(d)
                            if haskey(file, n*"_"*string(i))
                                delete_object(file, n*"_"*string(i))
                                file[n*"_"*string(i)] = d[i]
                            end
                        end
                    end
                end
            end            
        else
            try 
                h5open(path, "r+") do file
                    write(file, n, d)
                end
            catch
                if !isfile(path)
                    h5open(path, "w") do file
                        write(file, n, d) 
                    end
                else
                    h5open(path, "r+") do file
                        if haskey(file, n)
                            delete_object(file, n)
                            file[n] = d
                        end
                    end
                end
            end
        end
    end
end

function simple_save(struc, path::String)
    function Name(arg)
        string(arg)
    end
    fnames = fieldnames(typeof(struc))
    for fn in fnames 
        n = Name(fn)
        d = getfield(struc, fn)

        try 
            h5open(path, "r+") do file
                write(file, n, d)
            end
        catch
            h5open(path, "w") do file
                write(file, n, d) 
            end
        end
    end
end

parse_numbers(s, delimiter) = parse(Float64, split(s, delimiter, keepempty=false)[end])

function read_ITensor_vector(d, f, prefix, delimiter)
    # Read in the gates 
    Ts = []
    idx = []
    for key in keys(d)
        if startswith(key, prefix)
            T = read(f, key, ITensor)
            push!(Ts, T)
            push!(idx, parse_numbers(key, delimiter))
        end
    end
    return Ts[sortperm(idx)]
end

function load_params(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return Parameters(d["N"], d["Nx"], d["Ny"], d["Nsites"], 
            d["yperiodic"], d["doping"],d["max_phonons"], d["init_phonons"],
            d["μ"], d["εd"],d["εp"], d["tpd"], d["tpp"], d["Upd"],d["Upp"], d["Udd"],
            d["ω"], d["g0pp"],d["g0dd"], d["g1pd"], d["g1dp"], d["g1pp"],
            d["DMRG_numsweeps"], d["DMRG_noise"],d["DMRG_maxdim"],
            d["DMRG_cutoff"], d["DMRG_LBO"], d["max_LBO_dim"],d["min_LBO_dim"],
            d["mid"], d["T"], d["τ"], 
            d["TEBD_cutoff"], d["TEBD_maxdim"], d["TEBD_LBO"])
end

function load_equilibrium_correlations(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return EquilibriumCorrelations(d["start"], d["stop"], d["spin"], d["charge"], 
                                    d["sSC"], d["pSC"], d["dSC"])
end

function load_tebd_correlations(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return TEBDResults(d["entropy"], d["self_overlap"],
                                d["corrs"], nothing, nothing)
end

function load_dmrg_results_minimal(loadpath::String)
    f = h5open(loadpath,"r")
    d = read(f)
    return DMRGResults(nothing, d["ground_state_energy"], 
                            d["ground_state_entropy"], nothing,
                            d["charge_density"], d["phonon_density"], d["spin_density"])
end

function load_structs(loadpath::String; load_gs=true, 
                        load_phi=true, load_psi=true)
    f = h5open(loadpath,"r")
    d = read(f)

    # Read in the states 
    if load_gs
        ground_state = read(f, "ground_state", MPS)
        sites = siteinds(ground_state)
        hm = ThreeBandModel(params, sites)
    else
        ground_state = nothing
        sites = nothing
        hm = nothing
    end
    if load_phi
        phi_t = read(f, "phi_t", MPS)
    else
        phi_t = nothing
    end

    if load_psi
        psi_t = read(f, "psi_t", MPS)
    else
        psi_t = nothing
    end
 
    if d["DMRG_LBO"]
        # Read in the ITensor vectors
        optimized_basis = read_ITensor_vector(d, f, "optimized_basis_", "_")
    else
        optimized_basis = nothing
    end

    # Make sure to close the file after we are done reading 
    close(f)

    # Prepare the structs
    params = Parameters(d["N"], d["Nx"], d["Ny"], d["Nsites"], 
                        d["yperiodic"], d["doping"],d["max_phonons"], d["init_phonons"],
                        d["μ"], d["εd"],d["εp"], d["tpd"], d["tpp"], d["Upd"],d["Upp"], d["Udd"],
                        d["ω"], d["g0pp"],d["g0dd"], d["g1pd"], d["g1dp"], d["g1pp"],
                        d["DMRG_numsweeps"], d["DMRG_noise"],d["DMRG_maxdim"],
                        d["DMRG_cutoff"], d["DMRG_LBO"], d["max_LBO_dim"],d["min_LBO_dim"],
                        d["mid"], d["T"], d["τ"], 
                        d["TEBD_cutoff"], d["TEBD_maxdim"], d["TEBD_LBO"])
    dmrg_results = DMRGResults(ground_state, d["ground_state_energy"], 
                            d["ground_state_entropy"], optimized_basis,
                            d["charge_density"], d["phonon_density"], d["spin_density"])
    equilibrium_corr = EquilibriumCorrelations(d["start"], d["stop"], 
                                                d["spin"], d["charge"], 
                                                d["sSC"], d["pSC"], d["dSC"])
    # tebd_results = TEBDResults(d["entropy"], d["self_overlap"],
    #                             d["corrs"], phi_t, psi_t)

    return params, hm, dmrg_results, equilibrium_corr#, tebd_results
end