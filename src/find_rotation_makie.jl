
@info "loading packages"
using GadgetIO, GadgetUnits
using Printf
using ProgressMeter
using CairoMakie
using ColorSchemes
using SPHtoGrid
@info "done"


snap_file = "/gpfs/work/pn68va/di67meg/PaperRepos/WrongWayRelics/data/snap_057"

syn = read_block(snap_file, "SYNE", parttype=0)
sel = findall(syn .> -46.0)
syn = syn[sel]
sorted = sortperm(syn)
syn = syn[sorted]
pos = read_block(snap_file, "POS", parttype=0)[:, sel[sorted]]
center = [3261.2192, -2981.6863, 1397.6654]

pos_default = pos .- center

maximum(syn)

function plot_relic(pos_default, syn, rot1, rot2, rot3, plot_name)

    f = Figure(resolution=(4*800, 800))

    println("plot 1")
    ax1 = Axis(f[1, 1], aspect=1)
    scatter!(ax1, pos_default[2, :], pos_default[3, :], color=syn, marker='.', 
            colormap="magma", colorrange=(-46,-39))

    println("plot 2")
    ax2 = Axis(f[1, 2], aspect=1)
    pos = rotate_3D(pos_default, rot1...)
    scatter!(ax2, pos[1, :], pos[2, :], marker='.',
        color=syn, colormap="magma", colorrange=(-46, -39))

    println("plot 3")
    ax3 = Axis(f[1, 3], aspect=1)
    pos = rotate_3D(pos_default, rot2...)
    scatter!(ax3, pos[1, :], pos[2, :], marker='.',
        color=syn, colormap="magma", colorrange=(-46, -39))

    println("plot 4")
    ax4 = Axis(f[1, 4], aspect=1)
    pos = rotate_3D(pos_default, rot3...)
    scatter!(ax4, pos[1, :], pos[2, :], marker='.',
        color=syn, colormap="magma", colorrange=(-46, -39))


    println("saving")
    save(plot_name, f)

    nothing
end


rot1 = [105.0, 180.0, 90.0]
#             l/r
rot2 = [120.0, 180.0, 90.0]

rot3 = [150.0, 180.0, 90.0]

plot_name = "/gpfs/work/pn68va/di67meg/PaperRepos/WrongWayRelics/Plots/rotation2.png"
plot_relic(pos_default, syn, rot1, rot2, rot3, plot_name)

# pretty good
rot1 = [120.0, 180.0, 90.0]
rot1 = [0.0, 45.0, 180.0]
rot1 = [30.0, 45.0, 180.0]
rot1 = [45.0, 45.0, 180.0]

