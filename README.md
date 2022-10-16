# LMPTools.jl

[![DOI](https://zenodo.org/badge/343526561.svg)](https://zenodo.org/badge/latestdoi/343526561)

**Convenience and utility functions for [LongwaveModePropagator.jl](https://github.com/fgasdia/LongwaveModePropagator.jl) and other propagation models for the Earth-ionosphere waveguide.**

A brief summary of the exported functions are below, but use Julia's built-in help feature `?` for more information on using each of these functions.

## Transmitters

A dictionary of some VLF transmitters as `Transmitter` types can be obtained:

```julia
TRANSMITTER
```

Note that only each transmitter latitude, longitude, and frequency is correct. The power for each is 1000 W. 

---

A convenience function `range` is defined to compute the great circle distance on the WGS84 ellipsoid between `Transmitter` and `Receiver` types:

```julia
tx = TRANSMITTER[:NAA]
rx = Receiver("Boulder", 40.0, -105.3, 0.0, VerticalDipole())
greatcircledist = range(tx, rx)  # should return 3.146860697592529e6 meters
```

---

Although it's not exported, this package also extends `GeographicLib.GeodesicLine` for arguments `(tx::Transmitter, rx::Receiver)`.

## Ground

[LongwaveModePropagator.jl](https://github.com/fgasdia/LongwaveModePropagator.jl) exports `GROUND`, which is a dictionary of standard combinations of conductivity and relative permittivity.

This package includes a global ground conductivity map and functions which can be used to obtain the ground conductivity, relative permittivity, or `GROUND` index code. The map is the same one used by [LWPC](https://apps.dtic.mil/sti/citations/ADA350375), which is based on [(Morgan, 1986)](https://apps.dtic.mil/sti/citations/AD0675771). 

```julia
get_groundcode(lat, lon)
get_sigma(lat, lon)
get_epsilon(lat, lon)
```

---

There is also a function

```julia
groundsegments(tx::Transmitter, rx::Receiver; resolution=20e3, require_permittivity_change=false)
```

which can be used to generate a `Vector` of `Ground`s and a `Vector` of distances at which the ground changes between `tx` and `rx`. The keyword argument `resolution` is the step size in meters taken between `tx` and `rx` along the great circle path used to sample the ground conductivity map. The keyword argument `require_permittivity_change` can be set to `true` to require not only the conductivity but the relative permittivity to change in order to begin a new ground segment. This helps reduce the number of segments when changes in the ground are small.

## Magnetic field

This package provides convenience functions for interacting with the IGRF magnetic field model through `igrf13syn` from [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl) and the CHAOS-7 magnetic field through the Python package [ChaosMagPy](https://github.com/ancklo/ChaosMagPy).

```julia
igrf(tx::Transmitter, rx::Receiver, year, dists; alt=60e3)
```

returns a `Vector{BField}` along the great circle path from `tx` to `rx` at each distance of `dists` in meters.

It is more efficient to use the form

```julia
igrf(geoaz, lat, lon, year; alt=60e3)
```

if the magnetic field is required at a single location `lat`, `lon`, if the geodetic azimuth `geoaz` of the path from transmitter to receiver is already known.

The `chaos` function uses the complete CHAOS model (internal and external sources).

```julia
chaos(geoaz, lat, lon, year; alt=60e3)
chaos(geoaz, lats, lons, year; alt=60e3)  # special form for multiple lats, lons
chaos(tx, rx, year, dists; alt=60e3)
```

By default this package uses CHAOS-7.11. Alternate model coefficient MATLAB `.mat` files can be loaded using the function `load_CHAOS_matfile`. Subsequent calls to `chaos` will use the most recently loaded model coefficients.

## Ionospheres

Although not directly the ionosphere, the solar zenith angle is frequently needed for basic ionosphere models. This package provides

```julia
zenithangle(lat, lon, dt::DateTime, Δτ=nothing)
zenithangle(lat, lon, yr::Int, mo::Int, dy::Int, hr::Real, Δτ=nothing)
```

to compute the solar zenith angle in degrees using Algorithm 5 from [Grena, 2012. doi: 10.1016/j.solener.2012.01.024](https://doi.org/10.1016/j.solener.2012.01.024).

---

The function `isday(sza)` returns `true` if `sza` is less than 98 degrees.

---

```julia
ferguson(lat, sza, dt::DateTime) 
ferguson(lat, sza, month_num::Integer)
```

returns `(h′, β)` from the model given by [Ferguson, 1980](https://apps.dtic.mil/sti/citations/ADA085399).

---

```julia
flatlinearterminator(sza, hp_day=74, hp_night=86, b_day=0.3, b_night=0.5; sza_min=90, sza_max=100)
```

returns `(h′, β)` using a simple model for the terminator that is homogeneous on the day side and homogeneous on the night side with a linear transition between `sza_min` and `sza_max` solar zenith angle.

---

```julia
fourierperturbation(sza, coeff; a=deg2rad(45)/2)
```

returns a "perturbation" based on a Fourier series which can be used to model more complicated h′, β exponential ionospheres.

## Example with LongwaveModePropagator.jl

Here is an example using several LMPTools functions with LongwaveModePropagator.jl to model the electric field amplitude and phase at a receiver.

- The propagation path is segmented using `groundsegments`.
- The ionosphere is defined in h′ and β using Fourier perturbations (`fourierperturbation`) to the Ferguson model (`ferguson`).
- Earth's magnetic field comes from the `igrf` model.

```julia
using Dates, LongwaveModePropagator
using LongwaveModePropagator: QE, ME
using LMPTools
using GeographicLib  # using Pkg; Pkg.add("https://github.com/anowacki/GeographicLib.jl")

dt = DateTime(2021, 2, 1)

# Propagation path from NAA in Maine to Boulder, Colorado.
tx = TRANSMITTER[:NAA]
rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())

# For efficiency, precompute the geographic azimuth between the transmitter and receiver.
geoaz = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

# and precompute a line between the transmitter and receiver.
line = GeodesicLine(tx, rx)

# Fourier perturbation coefficients
hcoeff = (1.345, 0.668, -0.177, -0.248, 0.096)
bcoeff = (0.01, 0.005, -0.002, 0.01, 0.015)

# Use `groundsegments` from LMPTools to divide the propagation path into segments
grounds, dists = groundsegments(tx, rx; resolution=20e3)

# Preallocate a vector of `HomogeneousWaveguide` to construct a `SegmentedWaveguide`
wvgs = Vector{HomogeneousWaveguide{Species}}(undef, length(dists))
for i in eachindex(dists)
    dist = dists[i]
    wpt = forward(line, dist)

    bfield = igrf(geoaz, wpt.lat, wpt.lon, year(dt))

    sza = zenithangle(wpt.lat, wpt.lon, dt)
    h, b = ferguson(wpt.lat, sza, dt)
    h += fourierperturbation(sza, hcoeff)
    b += fourierperturbation(sza, bcoeff)

    species = Species(QE, ME, z->waitprofile(z, h, b), electroncollisionfrequency)

    wvgs[i] = HomogeneousWaveguide(bfield, species, grounds[i], dist)
end
wvg = SegmentedWaveguide(wvgs)

gs = GroundSampler(range(tx, rx), Fields.Ez)
E, amplitude, phase = propagate(wvg, tx, gs)  # field at the receiver
```
