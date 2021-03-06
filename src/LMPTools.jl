"""
    LMPTools

Convenience and utility functions for LongwaveModePropagator.jl.
"""
module LMPTools

using Rotations, StaticArrays
using HDF5, Dates
using GeographicLib, SatelliteToolbox
using LongwaveModePropagator

project_path(parts...) = normpath(@__DIR__, "..", parts...)

const GROUND_DATA = h5read(project_path("data", "conductivity_data.h5"), "/")

const TRANSMITTER = Dict(
    :NAA => Transmitter("NAA", 44.646394, -67.281069, 24e3),
    :NAU => Transmitter("NAU", 18.39875, -67.177433, 40.75e3),
    :NML => Transmitter("NML", 46.366014, -98.335544, 25.2e3),
    :NLK => Transmitter("NLK", 48.203611, -121.917056, 24.8e3),
    :NPM => Transmitter("NPM", 21.420, -158.151, 21.4e3)
)

export TRANSMITTER
export range
export get_ground, get_groundcode, get_epsilon, get_sigma, groundsegments
export igrf
export zenithangle, isday, ferguson, flatlinearterminator, smoothterminator, fourierperturbation


###
# Extend library functions

function GeographicLib.GeodesicLine(tx::Transmitter, rx::Receiver)
    return GeodesicLine(tx.longitude, tx.latitude; lon2=rx.longitude, lat2=rx.latitude)
end

"""
    range(tx::Transmitter, rx::Receiver)

Return the range (great-circle distance, but on the WGS-84 ellipsoid) in meters between
`tx` and `rx`.
"""
function Base.range(tx::Transmitter, rx::Receiver)
    return inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).dist
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

Return a `Vector{BField}` at each distance in `dists` in meters along the path from `tx`
to `rx` in fractional `year`.

If applying `igrf` to a single distance, it is recommended to use the `geoaz`, `lat`, `lon`
form of [`igrf`](@ref). 
"""
function igrf(tx::Transmitter, rx::Receiver, year, dists; alt=60e3)
    line = GeodesicLine(tx, rx)
    geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

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

Return the solar zenith angle in degrees for universal `dt` with `Δτ` difference in seconds
between terrestrial time and universal time. The position `lat`, `lon` should be in degrees.
"""
function zenithangle(lat, lon, dt::DateTime, Δτ=nothing)
    y, m, d = yearmonthday(dt)
    h = hour(dt) + minute(dt)/60 + second(dt)/3600 + millisecond(dt)/3600_000

    return zenithangle(lat, lon, y, m, d, h, Δτ)
end

"""
    isday(sza)

Return `true` if solar zenith angle `sza` is less than 98°, otherwise return `false`.
"""
function isday(sza)
    abs(sza) < 98 ? true : false
end

"""
    ferguson(lat, sza, month_num::Int) → (h′, β)
    ferguson(lat, sza, dt::DateTime)

Ferguson ionosphere for geographic latitude `lat` in degrees, solar zenith angle `sza` in
degrees, and month number `month_num` or `DateTime` `dt`.

If the month term of the ionosphere fit is not desired, specify `month_num = 0`.

The higher order terms in [^Ferguson1980] (sunspot number and magnetic absorption) are ignored.

# References

[Ferguson1980]: J. A. Ferguson, "Ionospheric profiles for predicting nighttime VLF/LF
    propagation," NOSC/TR-530, 1980.
"""
function ferguson(lat, sza, month_num::Integer)
    lat = deg2rad(lat)
    sza = deg2rad(sza)

    cos_lat = cos(lat)
    cos_sza = cos(sza)

    month_scalar = cospi(2*(month_num-0.5)/12)

    hprime = 74.37 - 8.097*cos_sza + 5.779*cos_lat - 1.213*month_scalar
    beta = 0.3849 - 0.1658*cos_sza - 0.08584*cos_lat + 0.1296*month_scalar

    return hprime, beta
end
ferguson(lat, sza, dt::DateTime) = ferguson(lat, sza, month(dt))

"""
    flatlinearterminator(sza, hp_day=74, hp_night=86, b_day=0.3, b_night=0.5) → (h′, β)

A simple terminator model with flat day and night ionospheres and a linear transition
between `sza_min` and `sza_max`. Solar zenith angle `sza` is in degrees.
"""
function flatlinearterminator(sza, hp_day=74, hp_night=86, b_day=0.3, b_night=0.5;
    sza_min=90, sza_max=100)

    if sza_min < sza < sza_max
        # Linear trend between day and night
        m_hprime = (hp_night - hp_day) / (sza_max - sza_min)
        m_beta = (b_night - b_day) / (sza_max - sza_min)

        hprime = m_hprime*(sza - sza_min) + hp_day
        beta = m_beta*(sza - sza_min) + b_day
    elseif sza < sza_min
        hprime = hp_day
        beta = b_day
    else  # sza > sza_max
        hprime = hp_night
        beta = b_night
    end

    return hprime, beta
end

"""
    flatlinearterminator(sza::AbstractArray; hp_day=74, hp_night=86, b_day=0.3, b_night=0.5) → (h′, β)

Given an array of `sza`, return the corresponding arrays of h′s and βs.
"""
function flatlinearterminator(sza::AbstractArray, hp_day=74, hp_night=86, b_day=0.3, b_night=0.5;
    sza_min=90, sza_max=100)

    # Linear trend between day and night
    m_hprime = (hp_night - hp_day) / (sza_max - sza_min)
    m_beta = (b_night - b_day) / (sza_max - sza_min)

    hprimes = similar(sza, Float64)
    betas = similar(hprimes)
    for i in eachindex(sza)
        if sza_min < sza[i] < sza_max  
            hprimes[i] = m_hprime*(sza[i] - sza_min) + hp_day
            betas[i] = m_beta*(sza[i] - sza_min) + b_day
        elseif sza[i] < sza_min
            hprimes[i] = hp_day
            betas[i] = b_day
        else  # sza > sza_max
            hprimes[i] = hp_night
            betas[i] = b_night
        end
    end

    return hprimes, betas
end

"""
    smoothterminator(sza, hp_day=74, hp_night=86, b_day=0.3, b_night=0.5;
        sza_min=90, sza_max=100, steepness=1) → (h′, β)
    
A smooth terminator model using a logistic curve from day to night ionospheres as a function
of solar zenith angle `sza`.

The `steepness` parameter `1` with the default values of `sza_min` and `sza_max`
approximately asymptotes to the day and night values by `sza_min` and `sza_max` and is
centered on their midpoint.
"""
function smoothterminator(sza, hp_day=74, hp_night=86, b_day=0.3, b_night=0.5;
    sza_min=90, sza_max=100, steepness=1)

    hp_night > hp_day || throw(ArgumentError("`hp_night` should be greater than `hp_day`."))
    b_night > b_day || throw(ArgumentError("`b_night` should be greater than `b_day`."))
    sza_max > sza_min || throw(ArgumentError("`sza_max` should be greater than `sza_min`."))

    hprime_range = hp_night - hp_day
    beta_range = b_night - b_day
    sza_mid = (sza_max - sza_min)/2 + sza_min

    # logistic curve
    hprime = hprime_range/(1 + exp(-steepness*(sza - sza_mid))) + hp_day
    beta = beta_range/(1 + exp(-steepness*(sza - sza_mid))) + b_day

    return hprime, beta
end

"""
    fourierperturbation(sza, coeff; a=deg2rad(45)/2)

Return the perturbation produced by the Fourier series:
```
pert = coeff[1] + coeff[2]*cos(sza*π/a) + coeff[3]*cos(sza*2π/a) + coeff[4]*sin(sza*π/a) + coeff[5]*sin(sza*2π/a)
```
for solar zenith angle `sza` in degrees.

This can be used to produce simulated ionospheres in conjunction with [`ferguson`](@ref).

# Examples
```julia-repl
julia> lat, lon = 45, -105
julia> dt = DateTime(2020, 1, 15);
julia> sza = zenithangle(lat, lon, dt);
julia> h, b = ferguson(lat, sza, dt);

julia> coeffH = (1.345, 0.668, -0.177, -0.248, 0.096);
julia> hpert = fourierperturbation(sza, coeffH)
1.7991484480844766

julia> h + hpert
79.5769471288202
```
"""
function fourierperturbation(sza, coeff; a=deg2rad(45)/2)
    sza = deg2rad(sza)

    pert = coeff[1] +
        coeff[2]*cospi(sza/a) +
        coeff[4]*sinpi(sza/a) +
        coeff[3]*cospi(sza*2/a) +
        coeff[5]*sinpi(sza*2/a)
    
    return pert
end

end  # module
