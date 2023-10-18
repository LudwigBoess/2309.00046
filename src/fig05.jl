include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
using GadgetIO, GadgetUnits
using Printf
using ProgressMeter
using DelimitedFiles
using Statistics
using Distributions
using StatsBase
using PyPlot, PyPlotUtility
using ColorSchemes
@info "done"


function plot_Inu_alpha_histograms(snap_range, z_range, plot_name)

    #colors = [1, 8, 15, 15]
    colors = ["purple", "mediumblue", "teal", "teal"]
    projections = ["xz", "xz", "xz", "yz"]
    linestyles = ["-", "-", "-", ":"]

    fig = get_figure(2.5)
    plot_styling!(axis_label_font_size=18)

    subplot(1,2,1)
    ax1 = gca()
    axis_ticks_styling!(ax1)

    ax1.set_yscale("log")
    ax1.set_xlim([-20, -16.0])
    ax1.set_ylim([1.e-5, 1.e-2])
    ax1.set_ylabel("Fraction of Image  " * L"N_\mathrm{pix}/N_\mathrm{tot}")
    ax1.set_xlabel("Synchrotron Intensity  " * L"\log_{10}(I_{ν,1.4 \mathrm{ GHz}})" * " [erg s" * L"^{-1}" * " Hz" * L"^{-1}" * " cm" * L"^{-2}" * "]")
    #ax1.xaxis.set_major_locator(plt.LinearLocator(10))

    subplot(1, 2, 2)
    ax2 = gca()
    axis_ticks_styling!(ax2)

    ax2.set_yscale("log")
    ax2.set_xlim([-1.5, -0.5])
    ax2.set_ylim([1.e-5, 1.e-2])
    ax2.set_ylabel("Fraction of Image  " * L"N_\mathrm{pix}/N_\mathrm{tot}")
    ax2.set_xlabel("Synchrotron Slope  " * L"\alpha_\mathrm{144 MHz}^\mathrm{1.4 GHz}")
    #ax2.xaxis.set_major_locator(plt.LinearLocator(10))


    for i = 1:length(snap_range)

        # c = ColorSchemes.batlow25[colors[i]]
        # c = (c.r, c.g, c.b)
        c = colors[i]

        data_file = data_path * "lower_relic_Inu_histograms_$(@sprintf("%03i", snap_range[i])).$(projections[i]).txt"
        data = readdlm(data_file)
        ax1.step(data[:, 1], data[:, 2] ./ 2048^2, where="post", linestyle=linestyles[i], color=c, label="z = $(@sprintf("%0.2f", z_range[i]))", lw=3)

        data_file = data_path * "lower_relic_alpha_histograms_$(@sprintf("%03i", snap_range[i])).$(projections[i]).txt"
        data = readdlm(data_file)
        ax2.step(data[:, 1], data[:, 2] ./ 2048^2, where="post", linestyle=linestyles[i], color=c, label="z = $(@sprintf("%0.2f", z_range[i])), $(projections[i])-plane", lw=3)
    end

    ax2.legend(frameon=false)

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end


snap_range = [56, 57, 58, 58]
z_range = [0.31515877622358834,
    0.29413466507090624,
    0.2734466450866466,
    0.2734466450866466]

plot_name = plot_path * "Fig05.pdf"

plot_Inu_alpha_histograms(snap_range, z_range, plot_name)