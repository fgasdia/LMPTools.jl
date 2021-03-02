"""
    LMPTools

Convenience and utility functions for LongwaveModePropagator.jl.
"""
module LMPTools

using Rotations, StaticArrays
using HDF5
using GeographicLib
using LongwaveModePropagator

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
export igrf


###
# Extend library functions

function GeographicLib.GeodesicLine(tx::Transmitter, rx::Receiver)
    return GeodesicLine(tx.longitude, tx.latitude; lon2=rx.longitude, lat2=rx.latitude)
end


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
    groundsegments(tx::Transmitter, rx::Receiver; resolution=20e3, require_permittivity_change=false)

Return a `Vector` of `Ground`'s and a `Vector` of distances at which the ground changes
between `tx` and `rx` using path `resolution` in meters.

If `require_permittivity_change == true`, then if the relative permittivity doesn't change
we don't consider it a ground change. This helps reduce the number of short small changes
that occur during a long path which have the same permittivity but a small change in
conductivity that has a small influence on the electric field in the EIWG but comes a the
cost of additional (potentially nearly identical) segments in the waveguide model.
"""
function groundsegments(tx::Transmitter, rx::Receiver; resolution=20e3, require_permittivity_change=false)

    line = GeodesicLine(tx, rx)
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

"""
    igrf(geoaz, lat, lon, year; alt=60e3)

Return a `BField` from IGRF-13 for position (`lat`, `lon`)°  at fractional
`year`. `geoaz`° is the geodetic azimuth of the path from transmitter to receiver.
By default, the magnetic field at an altitude of 60,000 meters is returned,
but this can be overridden with the `alt` keyword argument.
"""
function igrf(geoaz, lat, lon, year; alt=60e3)
    n, e, v, t = igrf13syn(0, year, 1, alt/1000, 90-lat, mod(lon, 360))  # in nT!, @60km altitude

    # Rotate the nev frame to the propagation path xyz frame
    # negate az to correct rotation direction for downward pointing v
    R = RotXZ(π, -deg2rad(geoaz))
    Rd = R*SVector(n,e,v)

    return BField(t*1e-9, Rd[1]/t, Rd[2]/t, Rd[3]/t)
end

"""
    igrf(tx::Transmitter, rx::Receiver, year, dists; alt=60e3)

Return a `Vector{BField}` at each distance in `dists` in meters along the path from position
(`lat1`, `lon1`)° and (`lat2`, `lon2`)° in fractional `year`.
"""
function igrf(tx::Transmitter, rx::Receiver, year, dists; alt=60e3)
    line = GeodesicLine(tx, rx)
    geoaz, _ = inverse(lon1, lat1, lon2, lat2)

    bfields = Vector{BField}(undef, length(dists))
    for (i, d) in enumerate(dists)
        lo, la, _ = forward(line, d)
        bfields[i] = igrf(geoaz, la, lo, year; alt=alt)
    end
    return bfields
end

"""
    igrf(tx::Transmitter, rx::Receiver, year, dists; alt=60e3)

Return a `Vector{BField}` at each distance in `dists` in meters along the path from `tx` to
`rx` in `year`.
"""
function igrf(tx::Transmitter, rx::Receiver, year, dists; alt=60e3)
    return igrf(tx.latitude, tx.longitude, rx.latitude, rx.longitude, year, dists; alt=alt)
end

end