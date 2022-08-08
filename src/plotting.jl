## IMPORTS ## 
using Plots
using FFTW
using NPZ
using Statistics
include(joinpath(@__DIR__,"model.jl"))

## EQUILIBRIUM CORRELATIONS ##

function plot_charge_density(dmrg_results::Union{DMRGResults,DMRGResultsMinimal})
    n = dmrg_results.charge_density'
    n = n[1:end-2,:] # remove the last rung
    plot(1:length(n[1:3:end,1]), n[1:3:end,:], label="py")
    plot!(1:length(n[2:3:end,1]), n[2:3:end,:], label="d")
    plot!(1:length(n[3:3:end,1]), n[3:3:end,:], label="px")
    ylabel!("⟨n⟩")
    xlabel!("Site")
    title!("Electron density")
end

function plot_phonon_density(dmrg_results::Union{DMRGResults,DMRGResultsMinimal}; ylims=nothing)
    n = dmrg_results.phonon_density'
    n = n[1:end-2,:] # remove the last rung
    plot(1:length(n[1:3:end,1]), n[1:3:end,:], label="py")
    plot!(1:length(n[2:3:end,1]), n[2:3:end,:], label="d")
    plot!(1:length(n[3:3:end,1]), n[3:3:end,:], label="px")
    if !isnothing(ylims)
        ylims!((ylim[0],ylim[1]))
    end
    ylabel!("⟨nb⟩")
    xlabel!("Site")
    title!("Phonon density")
end

function plot_spin_density(dmrg_results::Union{DMRGResults,DMRGResultsMinimal}; ylims=nothing)
    n = dmrg_results.spin_density'
    n = n[1:end-2,:] # remove the last rung
    plot(1:length(n[1:3:end,1]), n[1:3:end,:], label="py")
    plot!(1:length(n[2:3:end,1]), n[2:3:end,:], label="d")
    plot!(1:length(n[3:3:end,1]), n[3:3:end,:], label="px")
    if !isnothing(ylims)
        ylims!((ylim[0],ylim[1]))
    end
    ylabel!("⟨Sz⟩")
    xlabel!("Site")
    title!("Spin density")
end

function _plot_equilibrium_correlations(corrs, corrtype, start::Int, stop::Int)
    xrange = collect(1:stop-start+1)
    xrange = log10.(xrange)
    
    plot(xrange, corrs')
    title!(corrtype*"-"*corrtype*" correlation")
    xlabel!("Distance from centre site (log10)")
    ylabel!("Correlation (log10)")
end

function plot_equilibrium_correlations(eq_corrs::EquilibriumCorrelations, 
                                        corrtype::String)

    if corrtype=="spin"
        corrs = eq_corrs.spin
    elseif corrtype=="charge"
        corrs = eq_corrs.charge
    elseif corrtype=="sSC"
        corrs = eq_corrs.sSC
    elseif corrtype=="pSC"
        corrs = eq_corrs.pSC
    elseif corrtype=="dSC"
        corrs = eq_corrs.dSC  
    else
        @error "Invalid correlation type"
    end                               

    _plot_equilibrium_correlations(corrs, corrtype, eq_corrs.start, eq_corrs.stop)
end

## CHECKING TEBD RESULTS ##

function plot_entropy(tebd_results::TEBDResults)
    ent = tebd_results.entropy
    niters = size(ent)[2]
    ϕ_entropy = ent[1,:]
    ψ_entropy = ent[2,:]
    plot(1:niters, ϕ_entropy, label="ϕ(t)")
    plot!(1:niters, ψ_entropy, label="ψ(t)")
    title!("Von Neumann Entropy")
    xlabel!("Iteration")
end

function plot_overlap(tebd_results::TEBDResults)
    plot(1:length(tebd_results.self_overlap), tebd_results.self_overlap)
end

## TIME-DEPENDENT CORRELATIONS ## 

function plot_correlation_function(tebd_results::TEBDResults)
    heatmap(LinearAlgebra.norm.(tebd_results.corrs'),c=:heat)
    title!("Correlation function")
    xlabel!("Site")
    ylabel!("Time")
end

function make_spectral_fcn(corrs, p::Parameters; left_offset::Int=0)
    f1 = fftshift(fft(p.τ * corrs)')
    f2 = fftshift(fft(p.τ * reverse(corrs, dims=1))')

    qs = 2 * π * fftshift(fftfreq(p.N, 1))
    ωs = 2 * π * fftshift(fftfreq(size(corrs)[2], 1/p.τ))

    # Compensate for midpoint 
    sqw = zeros(size(f1))
    for i in 1:p.N 
        sqw[:,i] = imag.(exp.(1im * qs[i] * (p.mid-left_offset)) * f1[:,i])
                    + imag.(exp.(1im * qs[i] * (p.N-1-p.mid-left_offset)) * f2[:,i])
    end

    # If working with an even number of sites, must average across midline
    if iseven(p.N)
        sqw = (reverse(sqw, dims=2) + sqw) ./2
    end

    return real.(sqw)/π, ωs, qs 
end

function decay(u;shape="exponential")
    if shape=="exponential"
        t = collect(1:(length(u)))
        λ = 0.1/sqrt(length(t))
        v = exp.(-λ*t)
    else
        @error "Not implemented"
    end
    u .* v
end

function plot_spectral_function(tebd_results::TEBDResults, p::Parameters; 
                                smooth_signal=true, lims=nothing)
    corrs = tebd_results.corrs # num_time_steps x num_sites
    N, U, t, ω, g0, g1 = p.N, p.U, p.t, p.ω, p.g0, p.g1
    # Optionally convolve raw time data with decaying exponential?
    if smooth_signal
        # decaying exponential
        corrs = reverse(hcat(decay.(eachrow(corrs))...)',dims=2)
        # pad with zeros 
        corrs = hcat(corrs, zeros(size(corrs))) 
    end
    # Calculate the spectral function 
    ff, ωs, qs = make_spectral_fcn(corrs, p)

    # Zoom in to the relevant bit 
    if isnothing(lims)
        nstep = size(corrs)[2]
        lims = (floor(Int, nstep/2) - floor(Int, 0.07*nstep),floor(Int, nstep/2) + floor(Int, 0.07*nstep))
    end
    ff = ff[lims[1]:lims[2],:]
    ωs = ωs[lims[1]:lims[2]]

    # Plot 
    maxval = maximum(abs.(ff))
    heatmap(qs, ωs, abs.(ff), c=:bwr, clims=(-maxval, maxval))
    title!("N=$N, U=$U, t=$t, ω=$ω, g0=$g0, g1=$g1")
    xlabel!("Momentum")
    ylabel!("Frequency")
end

function plot_spectral_function_slice(tebd_results::TEBDResults, p::Parameters; slice=0)
    findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]
    
    ff, ωs, qs = make_spectral_fcn(tebd_results.corrs, p)
    ω = findnearest(ωs,slice)
    scatter(qs, abs.(ff[ω,:]),color="blue", label=nothing)
    plot!(qs, abs.(ff[ω,:]),color="blue", label=nothing)
    title!("Slice of S(q,ω) at ω=$slice")
end

function compare_to_ED(tebd_results::TEBDResults, p::Parameters)
    @assert p.τ == 0.01 # Time scale of the ED results
    corrs = tebd_results.corrs
    ED_corrs = npzread("/Users/nicole/Dropbox/Grad school/Devereaux lab/ITensors.jl/examples/dmrg/ed.npy")
    numsteps = size(corrs)[2]
    plot(1:numsteps, real.(corrs[4,:]), label="DMRG real part")
    plot!(1:numsteps, real.(ED_corrs[1:numsteps]), label="ED real part")
    plot!(1:numsteps, imag.(corrs[4,:]), label="DMRG complex part")
    plot!(1:numsteps, imag.(ED_corrs[1:numsteps]), label="ED complex part")
end