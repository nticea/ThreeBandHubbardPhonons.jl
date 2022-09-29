## IMPORTS ## 
using Plots
using FFTW
using NPZ
using Statistics
using CurveFit
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

function square_residual(y, ỹ)
    sum((ỹ .- y).^2)
end

function exponential_fit(x, y)
    a,b = exp_fit(x, y)
    fit_y = a.*exp.(b.*x)
    err = square_residual(fit_y, y)
    return a,b,fit_y,err
end

function power_law_fit(x, y)
    a,b = power_fit(x, y)
    fit_y = a.*(x.^b)
    err = square_residual(fit_y, y)
    return a,b,fit_y,err
end 

function compare_fits(x, y)
    exp_a, exp_b, exp_fit, exp_err  = exponential_fit(x,y)
    power_a, power_b, power_fit, power_err = power_law_fit(x,y)
    if power_err < exp_err
        return "power-law fit", power_a, power_b, power_fit
    else
        return "exponential fit", exp_a, exp_b, exp_fit
    end
end

function _plot_equilibrium_correlations(corrs, corrtype; do_fit=false, fit_type=nothing)
    corrs = corrs[1:end] # discard the correlation at zero distance (for now...)
    xrange = collect(1:length(corrs))
    #xrange = collect(0:length(corrs)-1)

    if do_fit
        if isnothing(fit_type)
            type, a, b, fit = compare_fits(xrange, abs.(corrs))
        else
            if fit_type=="exponential"
                a, b, fit, _ = exponential_fit(xrange, abs.(corrs))
                type = "exponential fit"
            elseif fit_type=="power"
                a, b, fit, _ = power_law_fit(xrange, abs.(corrs))
                type = "power-law fit"
            else
                @error "Fit type not recognized"
            end
        end
    end

    plot(log10.(xrange), log10.(abs.(corrs)))
    title!(corrtype*"-"*corrtype*" correlation")
    xlabel!("Distance from centre site (log10)")
    ylabel!("Correlation (log10)")
    plot!(log10.(xrange), log10.(fit), label=type)
end

function _plot_equilibrium_correlations_fitted(corrs, corrtype; fit_type=nothing)
    corrs = corrs[1:end] # discard the correlation at zero distance (for now...)
    xrange = collect(1:length(corrs))
    #xrange = collect(0:length(corrs)-1)

    if isnothing(fit_type)
        type, a, b, fit = compare_fits(xrange, abs.(corrs))
    else
        if fit_type=="exponential"
            a, b, fit, _ = exponential_fit(xrange, abs.(corrs))
            type = "exponential fit"
        elseif fit_type=="power"
            a, b, fit, _ = power_law_fit(xrange, abs.(corrs))
            type = "power-law fit"
        else
            @error "Fit type not recognized"
        end
    end

    plot(log10.(xrange), log10.(abs.(corrs)))
    title!(corrtype*"-"*corrtype*" correlation")
    xlabel!("Distance from centre site (log10)")
    ylabel!("Correlation (log10)")
    plot!(log10.(xrange), log10.(fit), label="$type, k=$b")
end

function _plot_equilibrium_correlations(corrs, corrtype)
    corrs = corrs[1:end] # discard the correlation at zero distance (for now...)
    xrange = collect(1:length(corrs))

    plot(log10.(xrange), log10.(abs.(corrs)))
    title!(corrtype*"-"*corrtype*" correlation")
    xlabel!("Distance from centre site (log10)")
    ylabel!("Correlation (log10)")
end

function plot_equilibrium_correlations(eq_corrs::EquilibriumCorrelations, 
                                        corrtype::String; do_fit=true, fit_type=nothing)

    if corrtype=="spin"
        corrs = eq_corrs.spin
    elseif corrtype=="charge"
        corrs = eq_corrs.charge
    elseif corrtype=="dSC_dxdx"
        corrs = eq_corrs.dSC_dxdx
    elseif corrtype=="dSC_dpx"
        corrs = eq_corrs.dSC_dpx
    elseif corrtype=="dSC_dydy"
        corrs = eq_corrs.dSC_dydy
    elseif corrtype=="dSC_pyd"
        corrs = eq_corrs.dSC_pyd
    elseif corrtype=="dSC_pypx"
        corrs = eq_corrs.dSC_pypx
    elseif corrtype=="dSC_py1px2"
        corrs = eq_corrs.dSC_py1px2
    else
        @error "Invalid correlation type"
    end

    if do_fit
        _plot_equilibrium_correlations_fitted(corrs, corrtype, fit_type=fit_type)
    else
        _plot_equilibrium_correlations(corrs, corrtype)
    end
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