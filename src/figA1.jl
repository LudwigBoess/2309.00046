include(joinpath(@__DIR__, "config.jl"))

using SpectralCRsUtility
using PyPlot, PyPlotUtility
using Printf
using ProgressMeter
#using ColorSchemes
using PyCall
cm = pyimport("cmasher")

"""
    read_cr_data(filename::String)

Reads the cooling spectrum from a binary file
"""
function read_cr_data(filename::String)

    f = open(filename, "r")
    Nfiles = read(f, Int64)
    Nbins = read(f, Int64)

    cr_norm = Vector{Vector{Float64}}(undef, Nfiles)
    cr_slope = Vector{Vector{Float64}}(undef, Nfiles)
    cr_cut = Vector{Float64}(undef, Nfiles)
    t = Vector{Float64}(undef, Nfiles)

    for i = 1:Nfiles
        cr_norm[i] = read!(f, Vector{Float64}(undef, Nbins))
        cr_slope[i] = read!(f, Vector{Float64}(undef, Nbins))
        cr_cut[i] = read(f, Float64)
        t[i] = read(f, Float64)
    end

    return cr_norm, cr_slope, cr_cut, t
end

"""
    get_synch_spectrum(cr_norm, cr_slope, cr_cut, ν_range)

Computes the synchrotron emissivity as a function of observational frequency.
"""
function get_synch_spectrum(cr_norm, cr_slope, cr_cut, ν_range)

    Lν = zeros(length(ν_range))

    Nbins = length(cr_norm)
    par = CRMomentumDistributionConfig(1.0, 1.e6, Nbins)

    for i = 1:length(ν_range)

        Lν[i] = synchrotron_emission(cr_norm, cr_slope, cr_cut, 5.e-6, par, ν0=ν_range[i],
            reduce_spectrum=true,
            integrate_pitch_angle=true,
            convert_to_mJy=false)
    end

    Lν
end

"""
    get_B_spectrum(cr_norm, cr_slope, cr_cut, B_range)

Computes the synchrotron emissivity as a function of magnetic field.
"""
function get_B_spectrum(cr_norm, cr_slope, cr_cut, B_range)

    Lν = zeros(length(B_range))

    Nbins = length(cr_norm)
    par = CRMomentumDistributionConfig(1.0, 1.e6, Nbins)

    for i = 1:length(B_range)

        Lν[i] = synchrotron_emission(cr_norm, cr_slope, cr_cut, B_range[i], par, ν0=1.4e9,
            reduce_spectrum=true,
            integrate_pitch_angle=true,
            convert_to_mJy=false)
    end

    Lν
end

"""
    α(q)

Analytic slope of a DSA injected synchrotron spectrum.
"""
α(q) = (q - 3) / 2

"""
    B0(q)

Analytic slope for scaling of synchrotron emission with magnetic field strength.
"""
B0(q) = α(q) + 1

"""
    plot_figA1(plot_name)
"""
function plot_figA1(plot_name)

    t_max = 203.20655731856036 # Myrs
    q0 = 4.5

    ν_range = 10.0 .^ LinRange(8, 11, 100)
    B_range = 10.0 .^ LinRange(-9, -3, 100)

    # read the data 
    filename = data_path * "cooling_spectra.dat"
    norm, slope, cut, t = read_cr_data(filename)

    sm = plt.cm.ScalarMappable(cmap=PyPlot.cm.viridis,
                 norm=plt.Normalize(vmin=0.0, vmax=t_max))
    sm.set_array([])

    fig = get_figure(3.5)
    plot_styling!()
    gs = plt.GridSpec(1, 4, hspace=0.0, wspace=0.3, width_ratios=[1, 1, 1, 0.05], figure=fig)

    subplot(get_gs(gs, 0, 0))
    ax = gca()
    ax.set_xlim([0.8e3, 1.2e6])
    ax.set_ylim([1.e-39, 1.e-25])
    ax.set_xscale("log")
    ax.set_yscale("log")
    axis_ticks_styling!(ax)

    xlabel("Dimensionless Momentum  " * L"\hat{p}" * " [" * L"(m_e \: c)^{-1}" * "]")
    ylabel("Distribution Function  " * L"f(\hat{p})" * " [arb. units]")

    @showprogress for i ∈ 1:101

        cr_bound, cr_norm = SpectralCRsUtility.norm_spectrum(norm[i], slope[i], cut[i], 1.0, 1.e6, 1.0, 4)

        ax.plot(cr_bound[1:end-1], cr_norm, c=sm.to_rgba(t[i]))

    end

    plot([1.e4, 1.e5], [1.e-29, 1.e-29 * 10.0^(-q0)], color="k", linestyle="--", linewidth=2)
    text(2.e4, 1.e-31, L"q_0 \: = " * "$(-q0)", rotation=-45)


    subplot(get_gs(gs, 0, 1))
    ax = gca()
    ax.set_xlim([1.e8, 1.e11])
    ax.set_ylim([1.e-29, 1.e-26])
    ax.set_xscale("log")
    ax.set_yscale("log")
    axis_ticks_styling!(ax)

    xlabel("Obs. Frequency  " * L"\nu" * " [Hz]")
    ylabel("Synch. Emissivity  " * L"j_\nu" * " [erg s" * L"^{-1}" * "Hz" * L"^{-1}" * "cm" * L"^{-3}" * "]")

    @showprogress for i ∈ 1:101

        Lν_ν = get_synch_spectrum(norm[i], slope[i], cut[i], ν_range)

        ax.plot(ν_range, Lν_ν, c=sm.to_rgba(t[i]))

    end

    plot([1.e9, 1.e10], [2.e-27, 2.e-27 * 10.0^(-α(q0))], color="k", linestyle="--", linewidth=2)
    text(2.e9, 1.e-27, L"\alpha_0 \: = " * "$(-α(q0))", rotation=-38)

    subplot(get_gs(gs, 0, 2))
    ax = gca()
    ax.set_xlim([1.e-9, 1.e-3])
    ax.set_ylim([1.e-34, 1.e-24])
    ax.set_xscale("log")
    ax.set_yscale("log")
    axis_ticks_styling!(ax)
    locmin = plt.LogLocator(base=10.0, subs=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), numticks=20)
    ax.xaxis.set_minor_locator(locmin)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    locmaj = matplotlib.ticker.LogLocator(base=10, numticks=12)
    ax.xaxis.set_major_locator(locmaj)

    xlabel("Magnetic Field  " * L"B" * " [G]")
    ylabel("Synch. Emissivity  " * L"j_\nu" * " [erg s" * L"^{-1}" * "Hz" * L"^{-1}" * "cm" * L"^{-3}" * "]")

    @showprogress for i ∈ 1:101

        Lν_B = get_B_spectrum(norm[i], slope[i], cut[i], B_range)

        ax.plot(B_range, Lν_B, c=sm.to_rgba(t[i]))

    end

    x = [5.e-8, 5.e-6]
    y = [1.e-30, 1.e-30 * 100^(B0(q0))]
    ax.plot()
    plot(x, y, color="k", linestyle="--", linewidth=2)
    text(4.e-8, 2.e-29, L"\alpha_0 + 1\: = " * "$(α(q0)+1)", rotation=45)

    subplot(get_gs(gs, 0, 3))
    cax = gca()
    cb = colorbar(sm, cax=cax, fraction=0.046)
    cb.set_label("Time  " * L"t" * " [Myr]")
    cb.ax.tick_params(
        direction="in",
        which="major",
        size=6, width=1
    )

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end

snap_range = 0:100
plot_name = plot_path * "FigA1.pdf"
plot_figA1(plot_name)



"""
    Error calculation
"""
ν_range = 10.0 .^ LinRange(8, 11, 100)
filename = data_path * "cooling_spectra.dat"
norm, slope, cut, t = read_cr_data(filename)
Lν_ν = get_synch_spectrum(norm[1], slope[1], cut[1], ν_range)


alpha_slope = (log10.(Lν_ν[1]) .- log10.(Lν_ν[end])) ./ (log10(ν_range[1]) - log10(ν_range[end]))

Δ = abs(alpha_slope + 0.75)/0.75 * 100

println("relative error in slope: Δα = ", Δ, " %")