include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using Printf
using ProgressMeter
using DelimitedFiles
using Statistics
using Distributions
using StatsBase
using PyPlot, PyPlotUtility
using ColorSchemes

"""
    get_histograms(bins, Q)

Computes the histograms for `M_s` of all shocked particles
"""
function get_histograms(bins, Q)

    # fit histograms
    hist = fit(Histogram, Q, bins)

    return hist.weights
end

"""
    read_mach_data(snap)

Reads all shocked particles.
"""
function read_mach_data(snap)

    snap_base = data_path * "snap_$(@sprintf("%03i", snap))"

    mach = read_block(snap_base, "MACH")

    sel = findall(mach .>= 1.0)

    return mach[sel]
end

"""
    construct_bins_and_centers(bin_min, bin_max, Nbins)

Returns the bins and their center points.
"""
function construct_bins_and_centers(bin_min, bin_max, Nbins)

    # construct bin boundaries
    bins = LinRange(bin_min, bin_max, Nbins + 1)

    # construct bin centers
    bin_centers = Vector{Float64}(undef, Nbins)

    for i = 1:Nbins
        bin_centers[i] = 0.5 * (bins[i] + bins[i+1])
    end

    bins, bin_centers
end

"""
    write_histograms(snap, Nbins=50)

Constructs the Mach number histograms of all shocked particles and writes them to a txt file.
"""
function write_histograms(snap, Nbins=50)

    # construct boundaries and centers for binning
    mach_bins, mach_centers = construct_bins_and_centers(1.0, 10.0, Nbins)

    # reads the simulation data
    mach = read_mach_data(snap)

    # construct histograms
    mach_binned = get_histograms(mach_bins, mach)

    # write into single array for saving
    data = [ mach_centers mach_binned ]

    # save as txt file
    data_file = data_path * "lower_relic_mach_histograms_$(@sprintf("%03i", snap)).txt"
    writedlm(data_file, data)
end

write_histograms(56)
write_histograms(57)
write_histograms(58)