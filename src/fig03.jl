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




function plot_mach_histograms(snap_range, z_range, plot_name)

#    colors = [1, 8, 15]
    colors = ["purple", "mediumblue", "teal"]


    fig = get_figure(1.0)
    plot_styling!()

    ax = gca()
    axis_ticks_styling!(ax)

    ax.set_yscale("log")
    ax.set_xlim([1.0, 10.0])
    ax.set_ylim([10, 2.e5])
    ax.set_ylabel("Shocked Particles")
    ax.set_xlabel("Sonic Mach Number " * L"\mathcal{M}_s")
    ax.xaxis.set_major_locator(plt.LinearLocator(10))

    ax.axvspan(1.0, 2.0, color="k", alpha=0.2)

    for i = 1:length(snap_range)

        data_file = data_path * "lower_relic_mach_histograms_$(@sprintf("%03i", snap_range[i])).txt"
        data = readdlm(data_file)

        # c = ColorSchemes.batlow25[colors[i]]
        # c = (c.r, c.g, c.b)
        c = colors[i]

        step(data[:, 1], data[:, 2], where="post", color=c, label="z = $(@sprintf("%0.2f", z_range[i]))", lw=3)
    end
    axvline(2.0, color="k", linestyle=":", alpha=0.3, label=L"\mathcal{M}_{s,\mathrm{crit.}}")

    legend(frameon=false)

    savefig(plot_name, bbox_inches="tight")
    close(fig)
end


snap_range = [56, 57, 58]
z_range = [0.31515877622358834,
    0.29413466507090624,
    0.2734466450866466]

plot_name = plot_path * "Fig03.pdf"

plot_mach_histograms(snap_range, z_range, plot_name)