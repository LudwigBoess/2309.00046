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

    files = [map_path * "g55_$(@sprintf("%03i", snap)).rot$rot_num.$filename.xy.fits"
             for rot_num ∈ collect(0:3), filename ∈ filenames]


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
    annotate_time = falses(Nrows * Ncols)
    time_labels = ["" for _ = 1:Nrows*Ncols]


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
    plot_name = plot_path * "FigA1.pdf"

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
