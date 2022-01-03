# LMPTools.jl

**Convenience and utility functions for [LongwaveModePropagator.jl](https://github.com/fgasdia/LongwaveModePropagator.jl) and other propagation models for the Earth-ionosphere waveguide.**

A brief summary of the exported functions are below, but use Julia's built-in help feature `?` for more information on using each of these functions.

## Transmitters

A dictionary of some VLF transmitters as `Transmitter` types can be obtained:

```julia
TRANSMITTER
```

Note that only each transmitter latitude, longitude, and frequency is correct. The power and altitude use the default values for LongwaveModePropagator's `Transmitter` type.

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

This package includes a global ground conductivity map (possibly the same one used by [LWPC](https://apps.dtic.mil/sti/citations/ADA350375), but admittedly the exact provenance of this file is unknown) and functions which can be used to obtain the ground conductivity, relative permittivity, or `GROUND` index code.

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

By default this package uses CHAOS-7.9. Alternate model coefficient MATLAB `.mat` files can be loaded using the function `load_CHAOS_matfile`. Subsequent calls to `chaos` will use the most recently loaded model coefficients.

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
