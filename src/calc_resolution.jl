include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using Printf
using PyPlotUtility
using Unitful, UnitfulAstro

# read parameters from used map
snap = 58
projection = "xz"
fi = map_path * "g55_$(@sprintf("%03i", snap)).synchSB_1.4GHz.$projection.fits"
map_synch, par_synch = PyPlotUtility.read_map_par(1, 1, [fi], nothing, nothing)

# read data from snapshot
fi = data_path * "snap_058"
h = read_header(fi)

# define cosmology
c = cosmology(h)

# calculate physical size of 1 arcsec at distance of redshift z=0.27
θ = 1.0u"arcsecond"
Δ = arcmin_to_kpc(c, θ, 0.27)

# calculate resulting resultion of pixels in arcsec
println("resolution: ", par_synch.pixelSideLength * 1.0u"kpc" / Δ, " ´´")