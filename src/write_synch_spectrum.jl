"""
    This file computes the synchrotron spectrum for all snapshots, which is quite expensive.
    It's shared-memory parallelized and I recommend running it with at least 16 cores in a screen session, or on in a slurm job.
"""


include(joinpath(@__DIR__, "config.jl"))

using GadgetIO, GadgetUnits
using ProgressMeter
using SpectralCRsUtility
using Base.Threads
using Printf

"""
    synch_power_per_frequency(data, GU, par, ν0)

Computes the total synchrotron luminosity in all particles for a given frequency
"""
function synch_power_per_frequency(data, GU, par, ν0)

    j_ν = Vector{Float64}(undef, size(data["CReN"], 2))

    @threads for i ∈ eachindex(data["MASS"])

        # convert to physical units
        norm = GU.CR_norm .* 10.0 .^ data["CReN"][:, i]
        slope = Float64.(data["CReS"][:, i])
        cut = Float64(data["CReC"][i])

        B = √(data["BFLD"][1, i]^2 +
              data["BFLD"][2, i]^2 +
              data["BFLD"][3, i]^2)

        # emissivity in [erg/s/Hz/cm^3] from SpectralCRsUtility.jl
        j_ν[i] = synchrotron_emission(norm, slope, cut, B, par; ν0)

        m = data["MASS"][i] * GU.m_cgs
        rho = data["RHO"][i] * GU.rho_cgs

        # convert to [W/Hz]
        j_ν[i] *= m / rho * 1.e-7

    end

    return sum(j_ν)

end


"""
    get_synch_power(data, h)

Computes the synchrotron spectrum over a number of frequency bins.
"""
function get_synch_power(data, h)

    @info "setting up spectra"

    GU = GadgetPhysical(h)

    # define observation frequency in Hz
    Nnu = 50
    ν_arr = 10.0 .^ (LinRange(7, 10, Nnu))
    P = Vector{Float64}(undef, Nnu)

    # cr setup 
    Nbins = size(data["CReN"], 1)
    par = CRMomentumDistributionConfig(0.1, 1.e5, Nbins)

    @info "calculating synch emissivity on $(nthreads()) threads"
    flush(stdout)
    flush(stderr)

    m = data["MASS"] .* GU.m_cgs
    rho = data["RHO"] .* GU.rho_cgs

    @showprogress for iNu = 1:Nnu

        # get total synch power per frequency
        P[iNu] = synch_power_per_frequency(data, GU, par, ν_arr[iNu])

        flush(stdout)
        flush(stderr)
    end # Nu

    @info "done"
    flush(stdout)
    flush(stderr)

    return ν_arr, P
end


function write_synch_spectrum(filename, ν_arr, P)

    f = open(filename, "w")

    write(f, Int64(length(ν_arr)))
    write(f, ν_arr)
    write(f, P)

    close(f)
end

function get_data(snap_base)

    blocks = ["POS", "HSML", "MASS", "RHO",
        "BFLD", "CReN", "CReS", "CReC"
    ]

    return Dict(block => read_block(snap_base, block, parttype=0) for block ∈ blocks )
end


function compute_synch_power(snap_range)


    for snap ∈ snap_range

        @info "reading data"
        flush(stdout)
        flush(stderr)

        snap_base = data_path * "snap_$(@sprintf("%03i", snap))"
        data = get_data(snap_base)
        h = read_header(snap_base)

        ν_arr, P = get_synch_power(data, h)

        P_filename = data_path * "lower_relic_spec_$(@sprintf("%03i", snap)).dat"
        write_synch_spectrum(P_filename, ν_arr, P)
    end

end

snap_range = 56:58
compute_synch_power(snap_range)


