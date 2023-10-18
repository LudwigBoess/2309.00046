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


function read_map_data(snap, projection)

    fi = map_path * "g55_$(@sprintf("%03i", snap)).synchSB_1.4GHz.$projection.fits"
    map_synch, par_synch = PyPlotUtility.read_map_par(1, 1, [fi], nothing, nothing)

    fi = map_path * "g55_$(@sprintf("%03i", snap)).synch_slope.$projection.fits"
    map_slope, par_slope = PyPlotUtility.read_map_par(1, 1, [fi], nothing, nothing)

    return log10.(vcat(map_synch...)), vcat(map_slope...)

end

"""
    write_histograms(snap, Nbins=50)

Constructs the Mach number histograms of all shocked particles and writes them to a txt file.
"""
function write_histograms(snap, projection, Nbins=50)

    # reads the map data
    map_synch, map_slope = read_map_data(snap, projection)

    # construct boundaries and centers for binning
    Inu_bins, Inu_centers = construct_bins_and_centers(-20.0, -16.0, Nbins)
    alpha_bins, alpha_centers = construct_bins_and_centers(-1.5, -0.5, Nbins)

    # construct histograms
    Inu_binned = get_histograms(Inu_bins, map_synch)

    # write into single array for saving
    data = [Inu_centers Inu_binned]
    # save as txt file
    data_file = data_path * "lower_relic_Inu_histograms_$(@sprintf("%03i", snap)).$projection.txt"
    writedlm(data_file, data)

    # construct histograms
    alpha_binned = get_histograms(alpha_bins, map_slope)

    # write into single array for saving
    data = [alpha_centers alpha_binned]
    # save as txt file
    data_file = data_path * "lower_relic_alpha_histograms_$(@sprintf("%03i", snap)).$projection.txt"
    writedlm(data_file, data)
end

write_histograms(56, "xz")
write_histograms(57, "xz")
write_histograms(58, "xz")
write_histograms(58, "yz")



