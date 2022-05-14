using Dates, Printf
using Plots, CSV
using LongwaveModePropagator
import LongwaveModePropagator as LMP
using GeographicLib
using LMPTools

include(joinpath("..", "test", "utils.jl"))

const GROUND_DATA = LMPTools.GROUND_DATA
const LAT = GROUND_DATA["lat"]
const LON = GROUND_DATA["lon"]

# Ground

tx = TRANSMITTER[:NAA]
rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())
distances = 0:100e3:range(tx, rx)
az = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude).azi

sigmas = [get_sigma(y, x) for y in LAT, x in LON]
heatmap(LON, LAT, sigmas, clims=(0, 0.11))

epsilons = [get_epsilon(y, x) for y in LAT, x in LON]
heatmap(LON, LAT, epsilons, clims=(0, 20))

pts = [forward(tx.longitude, tx.latitude, az, d) for d in distances]
heatmap(LON, LAT, sigmas, clims=(0, 0.11), xlims=(-110, -65), ylims=(35, 48))
scatter!(getfield.(pts, :lon), getfield.(pts, :lat))

# Magnetic field

# Compare to LWPC
tx = Transmitter{VerticalDipole}("fdtdnaa", 44.646, -67.281, VerticalDipole(), Frequency(24e3), 100e3)
rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())
grounds, distances = groundsegments(tx, rx; require_permittivity_change=false)
bfields = igrf(tx, rx, 2020, distances)

species = Species(LMP.QE, LMP.ME, z->waitprofile(z, 82, 0.55), electroncollisionfrequency)
ground = GROUND[10]
wvg = SegmentedWaveguide([HomogeneousWaveguide(bfields[i], species, ground, distances[i]) for i in 1:length(distances)])
sampler = GroundSampler(0:5e3:2000e3, Fields.Ez)
e, a, p = propagate(wvg, tx, sampler)

ld, la, lp = readlog(LMPTools.project_path("data", "bfield.log"))

plot(sampler.distance/1000, a, label="LMP",
     ylabel="Amplitude (dB)", xlabel="Range (km)", linewidth=1.5)
plot!(ld, la, label="LWPC", linewidth=1.5)

# Zenith angle

dt = DateTime(2021, 2, 1)
lats = 15:89
lons = -160:-60
szas = [zenithangle(la, lo, dt) for la in lats, lo in lons]
heatmap(lons, lats, szas,
        color=cgrad(:starrynight, [50, 70, 80, 85, 90, 95, 100, 110, 130]/180, rev=true), clims=(0, 180))


# Ionospheres

dt = DateTime(2021, 2, 1)
lats = 15:89
lons = -160:-60
x = [ferguson(la, zenithangle(la, lo, dt), dt) for la in lats, lo in lons]
heatmap(lons, lats, getindex.(x,1), color=:amp, clims=(68, 90))

hprimes, betas = flatlinearterminator(szas)
heatmap(lons, lats, hprimes, color=:amp, clims=(68, 90))
heatmap(lons, lats, betas, color=:tempo, clims=(0.2, 0.6))

ionos = smoothterminator.(szas)
heatmap(lons, lats, getindex.(ionos,1), color=:amp, clims=(68, 90))
heatmap(lons, lats, getindex.(ionos,2), color=:tempo, clims=(0.2, 0.6))

# Example

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
