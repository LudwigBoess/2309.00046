using SPHtoGrid
using Printf

include(joinpath(@__DIR__, "config.jl"))

"""
    get_slope_image(im_144, im_1440, vmin)

Computes the synchrotron slope per pixel between 144 MHz and 1.4 GHz.
`vmin` is an arbitrary cutoff to avoid computing pixels with spectra that are too steep to be relevant.
"""
function get_slope_image(im_144, im_1440, vmin)

    # safety cutoff
    im_144[im_144.<0.0] .= 0.0
    im_1440[im_1440.<0.0] .= 0.0

    im_144_cut = copy(im_144)
    im_144_cut[im_1440.<vmin] .= 0.0
    im_1440_cut = copy(im_1440)
    im_1440_cut[im_1440.<vmin] .= 0.0

    im_slope = (log10.(im_144_cut) .- log10.(im_1440_cut)) ./ (log10(144.e6) - log10(1.4e9))

    im_slope
end

"""
    write_slope_image(snap, orientation)


"""
function write_slope_image(snap, orientation)

    filenames = ["synchSB_144MHz", "synchSB_1.4GHz"]


    files = [map_path * "g55_$(@sprintf("%03i", snap)).$filename.$orientation.fits"
            for filename ∈ filenames]

    println("loading data")
    im_144, par, snap_num, units = read_fits_image(files[1])
    im_1440, par, snap_num, units = read_fits_image(files[2])

    println("constructing slope image")
    im_slope = get_slope_image(im_144, im_1440, 1.e-20)

    filename = map_path * "g55_$(@sprintf("%03i", snap)).synch_slope.$orientation.fits"
    println("saving slope image")
    write_fits_image(filename, im_slope,
                    par,
                    units = "[]",
                    snap = snap_num)

    println("done")

end

function write_slope_relics(snap, orientation)

    filenames = ["synchSB_144MHz", "synchSB_1.4GHz"]


    files = [map_path * "g55_$(@sprintf("%03i", snap)).$filename.$orientation.fits"
             for filename ∈ filenames]

    println("loading data")
    im_144, par, snap_num, units = read_fits_image(files[1])
    im_1440, par, snap_num, units = read_fits_image(files[2])

    println("constructing slope image")
    im_slope = get_slope_image(im_144, im_1440, 1.e-20)

    filename = map_path * "g55_$(@sprintf("%03i", snap)).synch_slope.$orientation.fits"
    println("saving slope image")
    write_fits_image(filename, im_slope,
        par,
        units="[]",
        snap=snap_num)

    println("done")

end

for snap = 56:58
    println("snap $snap")

    #write_slope_relics(snap, "xy")
    write_slope_relics(snap, "xz")
    write_slope_relics(snap, "yz")
end