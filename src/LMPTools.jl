"""
    LMPTools

Convenience and utility functions for LongwaveModePropagator.jl.
"""
module LMPTools

using Rotations, StaticArrays
using HDF5, Dates
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
export igrf, zenithangle


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


###
# Sun position

"""
    zenithangle(lat, lon, y, m, d, h, Δτ=nothing)

Return the solar zenith angle in degrees for year `y`, month `m`, day `d`, and hour `h` in
universal time (UT1). `Δτ` is the difference between terrestrial time and universal time.
The position in `lat`, `lon` is in degrees.

If `Δτ` is `nothing`, it will be computed by a linear extrapolation. Although this is not
accurate because `Δτ` is irregular, an error up to 30 seconds does not significantly affect
the result.

This function uses the calculation-optimized formulation of Algorithm 5 from

    Grena, R., “Five new algorithms for the computation of sun position from 2010 to 2110,”
    Solar Energy, 86, 2012, pp. 1323-1337. doi: 10.1016/j.solener.2012.01.024

which has a maximum error of 0.0027° (~10 arcsec) between year 2010 and 2110.

Refraction is not considered in this implementation.
"""
function zenithangle(lat, lon, yr::Int, mo::Int, dy::Int, hr::Real, Δτ=nothing)
    θ, ϕ = deg2rad(lon), deg2rad(lat)

    # Time scale computation
    if mo <= 2
        mo += 12
        yr -= 1
    end

    t = trunc(Int, 365.25*(yr - 2000)) + trunc(Int, 30.6001*(mo + 1)) -
        trunc(Int, 0.01*yr) + dy + 0.0416667*hr - 21958

    if isnothing(Δτ)
        Δτ = 96.4 + 0.00158*t
    end

    te = t + 1.1574e-5*Δτ

    # Heliocentric ecliptic longitude L
    ωa = 0.0172019715
    a = SVector(3.33024e-2, 3.512e-4, 5.2e-6)
    b = SVector(-2.0582e-3, -4.07e-5, -9e-7)
    s1, c1 = sincos(ωa*te)
    s2, c2 = 2*s1*c1, (c1 + s1)*(c1 - s1)
    s3, c3 = s2*c1 + c2*s1, c2*c1 - s2*s1
    s, c = SVector(s1, s2, s2), SVector(c1, c2, c3)

    β = 2.92e-5
    dβ = -8.23e-5
    ω = SVector(1.49e-3, 4.31e-3, 1.076e-2, 1.575e-2, 2.152e-2, 3.152e-2, 2.1277e-1)
    d = SVector(1.27e-5, 1.21e-5, 2.33e-5, 3.49e-5, 2.67e-5, 1.28e-5, 3.14e-5)
    φ = SVector(-2.337, 3.065, -1.533, -2.358, 0.074, 1.547, -0.488)

    L = 1.7527901 + 1.7202792159e-2*te
    for k = 1:3
        L += a[k]*s[k] + b[k]*c[k]
    end
    L += dβ*s1*sin(β*te)
    for i = 1:7
        L += d[i]*sin(ω[i]*te + φ[i])
    end

    # Nutation correction
    ωn = 9.282e-4
    ν = ωn*te - 0.8
    sν, cν = sincos(ν)
    Δλ = 8.34e-5*sν
    λ = L + π + Δλ
    ϵ = 4.089567e-1 - 6.19e-9*te + 4.46e-5*cν
    sλ, cλ = sincos(λ)
    sϵ, cϵ = sincos(ϵ)

    α = atan(sλ*cϵ, cλ)
    α < 0 && (α += 2π)  # [-π/2, π/2] → [0, 2π]

    δ = asin(sλ*sϵ)
    H = 1.7528311 + 6.300388099*t + θ - α + 0.92*Δλ
    H = mod2pi(H + π) - π  # hour angle

    # Calculate zenith and azimuth angle with parallax correction for topocentric coords
    sϕ, cϕ = sincos(ϕ)
    sδ, cδ = sincos(δ)
    cH = cos(H)
    # sH, cH = sincos(H)

    se0 = sϕ*sδ + cϕ*cδ*cH
    ep = asin(se0) - 4.26e-5*sqrt(1 - se0^2)

    # Γ = atan(sH, cH*sϕ - sδ*cϕ/cδ)  # azimuth

    # Refraction correction to elevation
    # Δre = 0.08422*P/(273 + T)/tan(ep + 0.003138/(ep + 0.08919))

    z = π/2 - ep  # - Δre  # zenith angle

    return rad2deg(z)
end

"""
    zenithangle(lat, lon, dt::DateTime, Δτ=nothing)

Return the solar zenith angle in radians for universal `dt` with `Δτ` difference in seconds
between terrestrial time and universal time. The position `lat`, `lon` should be in degrees.
"""
function zenithangle(lat, lon, dt::DateTime, Δτ=nothing)
    y, m, d = yearmonthday(dt)
    h = hour(dt) + minute(dt)/60 + second(dt)/3600 + millisecond(dt)/3600_000

    return zenithangle(lat, lon, y, m, d, h, Δτ)
end

end