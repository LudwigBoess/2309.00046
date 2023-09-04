include(joinpath(@__DIR__, "config.jl"))

@info "loading packages"
# load packages
using GadgetIO
using PyPlot, PyPlotUtility
using Printf
using PyCall
cm = pyimport("cmasher")
@info "done"

function plot_relic_evolution()

    @info "plotting..."

    filenames = ["Xray", "CReE", "synchSB_1.4GHz", "synch_slope"]

    Ncols = 4
    Nrows = 4

    files1 = [map_path * "g55_$(@sprintf("%03i", snap)).$filename.yz.fits"
              for snap ∈ [56, 57, 58], filename ∈ filenames]

    files2 = [map_path * "g55_$(@sprintf("%03i", 58)).$filename.xz.fits" for filename ∈ filenames]

    files = Matrix{String}(undef, 4, 4)

    files[1:3, :] = files1
    files[4, :] = files2

    # limits for colorbars
    vmin_arr = [1.e-5, 1.e6, 1.e-20, -1.5]
    vmax_arr = [5.e-1, 2.e9, 5.e-17, -0.5]

    # colormaps per column
    im_cmap = [cm.sunburst, cm.freeze, "cubehelix", "Spectral_r"]

    # labels for colorbars
    cb_labels = [L"SB_{x, 0.1-2.4 \mathrm{ keV}}" * "[erg s" * L"^{-1}" * " cm" * L"^{-2}" * "]",
                L"E_{\mathrm{CR}e^- > 1 \mathrm{ GeV}}" * " [erg cm" * L"^{-2}" * "]",
                L"I_{ν,1.4 \mathrm{ GHz}}" * " [erg s" * L"^{-1}" * " Hz" * L"^{-1}" * " cm" * L"^{-2}" * "]",
                L"\alpha_\mathrm{144 MHz}^\mathrm{1.4 GHz}"]

    
    # redshift annotation
    annotate_time = trues(Nrows * Ncols)
    time_labels = ["" for _ = 1:Nrows*Ncols]
    time_labels[1] = "z = $(@sprintf("%0.2f", 0.31515877622358834))"
    time_labels[2] = "z = $(@sprintf("%0.2f", 0.29413466507090624))"
    time_labels[3] = "z = $(@sprintf("%0.2f", 0.2734466450866466 ))"
    time_labels[4] = "z = $(@sprintf("%0.2f", 0.2734466450866466 ))"

    # logarithmic colorbar
    log_map = trues(Nrows * Ncols)
    log_map[4] = false
    log_map[8] = false
    log_map[12] = false
    log_map[16] = false

    # scale line
    annotate_scale = falses(Nrows * Ncols)
    annotate_scale[1] = true
    annotate_scale[2] = true
    annotate_scale[3] = true
    annotate_scale[4] = true

    # scale label
    scale_kpc = 500.0
    scale_label = "500 kpc"

    # name of the plot
    plot_name = plot_path * "Fig02.pdf"

    plot_image_grid(Nrows, Ncols, files, im_cmap, cb_labels,
        vmin_arr, vmax_arr, plot_name,
        cutoffs=vmin_arr,
        mask_bad=trues(Nrows * Ncols),
        upscale=0.9,
        cb_label_offset=0.6,
        dpi=200;
        time_labels, annotate_time,
        log_map,
        annotate_scale,
        scale_kpc,
        scale_label
    )

    GC.gc()
end


plot_relic_evolution()
