"""
    LMPTools

Convenience and utility functions for LongwaveModePropagator.jl.
"""
module LMPTools

using Rotations, StaticArrays
using HDF5
using GeographicLib
using LongwaveModePropagator
import LongwaveModePropagator: BField

project_path(parts...) = normpath(@__DIR__, "..", parts...)

include("igrf13syn.jl")
const GROUND_DATA = h5read(project_path("data", "conductivity_data.h5"), "/")

const TRANSMITTER = Dict(
    :NAA => Transmitter("NAA", 44.646394, -67.281069, 24e3),
    :NAU => Transmitter("NAU", 18.39875, -67.177433, 40.75e3),
    :NML => Transmitter("NML", 46.366014, -98.335544, 25.2e3),
    :NLK => Transmitter("NLK", 48.203611, -121.917056, 24.8e3),
    :NPM => Transmitter("NPM", 21.420, -158.151, 21.4e3)
)

export TRANSMITTER
export get_ground, get_groundcode, get_epsilon, get_sigma, groundsegments

###
# Ground

function searchsortednearest(a, x)
    @assert issorted(a) "`a` must be sorted!"

    idx = searchsortedfirst(a, x)

    idx == 1 && return idx
    idx > length(a) && return length(a)
    a[idx] == x && return idx
    abs(a[idx]-x) < abs(a[idx-1]-x) ? (return idx) : (return idx-1)
end

"""
    get_ground(param, lat, lon)

Return the ground `param` "codemap", "epsilonmap", or "sigmamap" at position `lat`, `lon`
in degrees.

Convenience functions are available for the most common fields with arguments `lat`, `lon`
only:
    
    * [`get_groundcode`](@ref)
    * [`get_sigma`](@ref)
    * [`get_epsilon`](@ref)
"""
function get_ground(param, lat, lon)
    latidx = searchsortednearest(GROUND_DATA["lat"], lat)
    lonidx = searchsortednearest(GROUND_DATA["lon"], lon)

    return GROUND_DATA[param][latidx, lonidx]
end
get_groundcode(lat, lon) = get_ground("codemap", lat, lon)
get_sigma(lat, lon) = get_ground("sigmamap", lat, lon)
get_epsilon(lat, lon) = get_ground("epsilonmap", lat, lon)


"""
    groundsegments(lat1, lon1, lat2, lon2; resolution=20e3, require_permittivity_change=false)

Return a `Vector` of `Ground`'s and a `Vector` of distances at which the ground changes
between `lat1, lon1` and `lat2, lon2` using path `resolution` in meters.

If `require_permittivity_change == true`, then if the relative permittivity doesn't change
we don't consider it a ground change. This helps reduce the number of short small changes
that occur during a long path which have the same permittivity but a small change in
conductivity that has a small influence on the electric field in the EIWG but comes a the
cost of additional (potentially nearly identical) segments in the waveguide model.
"""
function groundsegments(lat1, lon1, lat2, lon2;
    resolution=20e3, require_permittivity_change=false)

    line = GeodesicLine(lon1, lat1; lon2=lon2, lat2=lat2)
    pts = waypoints(line, dist=resolution)

    lat, lon, dist = first(pts).lat, first(pts).lon, first(pts).dist
    lastground = GROUND[get_groundcode(lat, lon)]

    # initialize
    grounds = [lastground]
    distances = [dist]  # == 0.0

    for pt in pts
        lat, lon, dist = pt.lat, pt.lon, pt.dist
        ground = GROUND[get_groundcode(lat, lon)]

        if require_permittivity_change
            if ground.ϵᵣ != lastground.ϵᵣ
                push!(grounds, ground)
                push!(distances, dist)
                lastground = ground
            end
        elseif ground != lastground
            push!(grounds, ground)
            push!(distances, dist)
            lastground = ground
        end
    end

    return grounds, distances
end


###
# Magnetic field

end