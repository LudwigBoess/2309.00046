println("allocating cores")
using Distributed, ClusterManagers

# automatically decide if it needs to be run in slurm envirnoment
try
    println("allocating $(ENV["SLURM_NTASKS"]) slurm tasks")
    addprocs_slurm(parse(Int64, ENV["SLURM_NTASKS"]))
catch err
    if isa(err, KeyError)
        println("allocating 4 normal tasks")
        addprocs(4)
    end
end

println("done")
flush(stdout);
flush(stderr);

include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using SPHKernels, SPHtoGrid
using Printf
using ProgressMeter
using Base.Threads
using SpectralCRsUtility
using DelimitedFiles

"""
    get_synchrotron(data, nu, GU)

Computes the synchrotron emissivity for each particle in units [erg/s/Hz/cm^3]
"""
function get_synchrotron(data, nu, GU)

    println("synchrotron running on $(nthreads()) threads")

    Npart = length(data["CReC"])
    j_ν = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)

    @threads for i ∈ eachindex(j_ν)

        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]

        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        B = sqrt(data["BFLD"][1, i]^2 + data["BFLD"][2, i]^2 + data["BFLD"][3, i]^2)

        j_ν[i] = synchrotron_emission(norm, slope, cut, B, par, ν0=nu)
    end

    j_ν
end

"""
    get_CReE(data, GU)

Computes the total energy contained in electrons above 1GeV in [erg].
"""
function get_CReE(data, GU)

    println("CReE running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)


    Npart = length(data["CReC"])
    CReE = Vector{Float64}(undef, Npart)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)
    bounds = momentum_bin_boundaries(par)

    @threads for i = 1:Npart

        norm = 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        CReE[i] = cr_energy_in_range(norm, slope, cut, 1.0, bounds)
    end

    # convert to erg/cm^3
    return CReE .* (GU.E_cgs / GU.x_cgs^3)
end


"""
    make_quantity_maps(snap)

Writes all maps that are required for making Fig. 2.
"""
function make_quantity_maps(snap)

    println("running on $(nthreads()) threads")
    flush(stdout)
    flush(stderr)

    println("reading data")
    flush(stdout)
    flush(stderr)

    snap_base = data_path * "snap_$(@sprintf("%03i", snap))"

    h = read_header(snap_base)

    GU = GadgetPhysical(h)

    blocks = ["POS", "HSML", "MASS", "RHO", "U",
        "BFLD",
        "CReN", "CReS", "CReC"
    ]
    data = Dict(block => read_block(snap_base, block, parttype=0) for block ∈ blocks)

    println("done")
    flush(stdout)
    flush(stderr)

    image_path = map_path * "g55_$(@sprintf("%03i", snap))."

    if snap == 56
        center = [3.5, -3.2, 1.5] .* 1.e3 ./ GU.x_physical
        width = 1000.0 ./ GU.x_physical
    elseif snap == 57
        center = [3.5, -3.2, 1.5] .* 1.e3 ./ GU.x_physical
        width = 1000.0 ./ GU.x_physical
    elseif snap == 58
        center = [3.39, -2.7506, 1.467] .* 1.e3
        width = 1000.0 ./ GU.x_physical
    end


    # select kernel
    kernel = WendlandC4(Float64, 2)

    # convert to physical code units for mapping
    pos = data["POS"] .* GU.x_physical
    hsml = data["HSML"] .* GU.x_physical
    rho = data["RHO"] .* GU.rho_physical
    mass = data["MASS"] .* GU.m_physical

    xy_size = width
    z_size = width

    # define mapping parameters
    param = mappingParameters(center=center .* GU.x_physical,
        x_size=xy_size * GU.x_physical,
        y_size=xy_size * GU.x_physical,
        z_size=z_size * GU.x_physical,
        Npixels=1024)


    println("rho")
    flush(stdout); flush(stderr);
    rho_cgs = data["RHO"] .* GU.rho_cgs
    weights = part_weight_physical(length(rho_cgs), param)

    image_prefix = image_path * "rho"
    # weights = part_weight_physical(length(rho_cgs), param)
    # map_it(pos, hsml, mass, rho, rho_cgs, weights, 
    #         units="g/cm^2", reduce_image=false; 
    #         kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, rho_cgs, weights, 
            units="g/cm^2", reduce_image=false,
            projection="xz"; 
            kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, rho_cgs, weights, 
            units="g/cm^2", reduce_image=false,
            projection="yz";  
            kernel, snap, param, image_prefix)


    println("Xray")
    flush(stdout); flush(stderr);

    image_prefix = image_path * "Xray"
    T_keV = data["U"] .* GU.T_eV .* 1.e-3
    rho_cgs = data["RHO"] .* GU.rho_cgs
    Xray = x_ray_emissivity(T_keV, rho_cgs)

    weights = part_weight_physical(length(rho_cgs), param)

    # map_it(pos, hsml, mass, rho, Xray, weights, 
    #         units="erg/s/cm^2", reduce_image=false; 
    #         kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, Xray, weights, 
            units="erg/s/cm^2", reduce_image=false,
            projection="xz"; 
            kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, Xray, weights, 
            units="erg/s/cm^2", reduce_image=false,
            projection="yz";  
            kernel, snap, param, image_prefix)


    println("CReE")
    flush(stdout); flush(stderr);

    image_prefix = image_path * "CReE"
    CReE = get_CReE(data, GU)
    weights = part_weight_physical(length(rho_cgs), param)

    # map_it(pos, hsml, mass, rho, CReE, weights,
    #     units="erg/cm^2",
    #     reduce_image=false;
    #     kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, CReE, weights,
            units="erg/cm^2",
            reduce_image=false,
            projection="xz"; 
            kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, CReE, weights,
        units="erg/cm^2",
        reduce_image=false,
        projection="yz";
        kernel, snap, param, image_prefix)

    CReE = nothing
    GC.gc()


    println("synch SB")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "synchSB_144MHz"
    j_nu = get_synchrotron(data, 144.0e6, GU)

    weights = part_weight_physical(length(rho_cgs), param)

    # map_it(pos, hsml, mass, rho, j_nu, weights,
    #     units="erg/s/Hz/cm^2",
    #     reduce_image=false,
    #     projection="xy";
    #     kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, j_nu, weights,
        units="erg/s/Hz/cm^2",
        reduce_image=false,
        projection="xz";
        kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, j_nu, weights,
        units="erg/s/Hz/cm^2",
        reduce_image=false,
        projection="yz";
        kernel, snap, param, image_prefix)
    j_nu = nothing

    println("synch SB")
    flush(stdout)
    flush(stderr)

    image_prefix = image_path * "synchSB_1.4GHz"
    j_nu = get_synchrotron(data, 1.4e9, GU)
    weights = part_weight_physical(length(rho_cgs), param)

    # map_it(pos, hsml, mass, rho, j_nu, weights,
    #     units="erg/s/Hz/cm^2",
    #     reduce_image=false,
    #     projection="xy";
    #     kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, j_nu, weights,
        units="erg/s/Hz/cm^2",
        reduce_image=false,
        projection="xz";
        kernel, snap, param, image_prefix)
    map_it(pos, hsml, mass, rho, j_nu, weights,
        units="erg/s/Hz/cm^2",
        reduce_image=false,
        projection="yz";
        kernel, snap, param, image_prefix)
    j_nu = nothing

    data = nothing 
    GC.gc()

end


for snap = 56:58
    make_quantity_maps(snap)
end

