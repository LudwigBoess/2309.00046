include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using PyPlot, PyPlotUtility
using Printf
using SPHtoGrid
using Statistics
using ProgressMeter
using ColorSchemes


"""
    read_spectrum_sim(filename)

Reads the custom binary format for the relic spectrum.
"""
function read_spectrum_sim(filename)

    f = open(filename, "r")
    Nnu = read(f, Int64)

    ν_arr = read!(f, Vector{Float64}(undef, Nnu))
    P_sim = read!(f, Vector{Float64}(undef, Nnu))

    close(f)

    return ν_arr, P_sim
end

"""
    plot_synch_spectra_evolution(plot_name)

Plots Figure 4.
"""
function plot_synch_spectra_evolution(plot_name)

    # relevant snapshots
    snaps = [56, 57, 58]

    # redshifts of snapshots 
    z = [0.31515877622358834, 0.29413466507090624, 0.2734466450866466]

    # colors in colorscale array
    colors = [1, 8, 12]

    # helper function to get my plot styling
    fig = get_figure(1.0)
    plot_styling!()

    ax = gca()
    axis_ticks_styling!(ax)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim([10.0, 1.e4])
    ax.set_ylim([1.e22, 3.e25])
    ax.set_xlabel("Obs. Frequency  " * L"\nu" * " [MHz]")
    ax.set_ylabel("Synchrotron Luminosity  " * L"P(\nu)" * " [W Hz" * L"^{-1}" * "]")

    # plot spectral fits "by hand"
    # z = 0.27
    P0 = [8e23, 5.e22]
    ν_example = [1.3299518105153717e1, 7.196856730011514e1, 3.906939937054621e2, 1.2067926406393263e3]
    α = [0.735, 1.186]
    plot(ν_example[1:2], [P0[1], P0[1] * (ν_example[1] / ν_example[2])^α[1]], color="k", linestyle="--")
    plot(ν_example[3:4], [P0[2], P0[2] * (ν_example[3] / ν_example[4])^α[2]], color="k", linestyle="--")
    text(1.5e1, 2e23, L"\alpha_1 = -0.735", rotation=329, fontsize=12)
    text(3.5e2, 1.3e22, L"\alpha_2 = -1.18", rotation=315, fontsize=12)

    # z = 0.29
    P0 = 7e23
    ν_example = [1.e2, 1.e3]
    α = [0.66]
    plot(ν_example, [P0, P0 * (ν_example[1] / ν_example[2])^α[1]],
        color="k", linestyle="--")
    text(1.5e2, 1.9e23, L"\alpha = -0.66", rotation=330, fontsize=12)

    # z = 0.32
    P0 = [1.1e25, 4.e23]
    ν_example = [1.3299518105153717e1, 7.196856730011514e1, 1.2067926406393263e3, 6.551285568595523e3]
    α = [0.677, 1.0]
    plot(ν_example[1:2], [P0[1], P0[1] * (ν_example[1] / ν_example[2])^α[1]], color="k", linestyle="--")
    plot(ν_example[3:4], [P0[2], P0[2] * (ν_example[3] / ν_example[4])^α[2]], color="k", linestyle="--")
    text(2.0e1, 5.e24, L"\alpha_1 = -0.677", rotation=330, fontsize=12)
    text(2.0e3, 1.5e23, L"\alpha_2 = -1.0", rotation=320, fontsize=12)


    # read and plot all spectra
    for i = 1:3
        # spectrum in particles 
        filename = data_path * "lower_relic_spec_$(@sprintf("%03i", snaps[i])).dat"
        ν_part, P_part = read_spectrum_sim(filename)
        c = ColorSchemes.batlow25[colors[i]]
        plot(ν_part .* 1.e-6, P_part, label="z = $(@sprintf("%0.2f", z[i]))", lw=3, color=(c.r, c.g, c.b))
    end

    legend(frameon=false)

    savefig(plot_name, bbox_inches="tight", dpi=200)
    close(fig)

end

plot_name = plot_path * "Fig04.pdf"
plot_synch_spectra_evolution(plot_name)

