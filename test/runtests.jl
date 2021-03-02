using Test, LMPTools
using CSV
using GeographicLib
using LongwaveModePropagator
const LMP = LongwaveModePropagator

include("lwpc_utils.jl")

const GROUND_DATA = LMPTools.GROUND_DATA
const LAT = GROUND_DATA["lat"]
const LON = GROUND_DATA["lon"]


@testset "LMPTools" begin
    # Ground
    @test extrema(GROUND_DATA["lat"]) == (-90, 90)
    @test extrema(GROUND_DATA["lon"]) == (-180, 180)

    @test LMPTools.searchsortednearest(GROUND_DATA["lat"], 78) ==
        findfirst(isequal(78), GROUND_DATA["lat"])
    @test LMPTools.searchsortednearest(GROUND_DATA["lat"], 78.1) ==
        findfirst(isequal(78), GROUND_DATA["lat"])
    @test LMPTools.searchsortednearest(GROUND_DATA["lat"], 78.4) ==
        findfirst(isequal(78.5), GROUND_DATA["lat"])
    @test LMPTools.searchsortednearest(GROUND_DATA["lat"], -78) ==
        findfirst(isequal(-78), GROUND_DATA["lat"])
    @test LMPTools.searchsortednearest(GROUND_DATA["lat"], -78.1) ==
        findfirst(isequal(-78), GROUND_DATA["lat"])
    @test LMPTools.searchsortednearest(GROUND_DATA["lat"], -78.4) ==
        findfirst(isequal(-78.5), GROUND_DATA["lat"])

    for la in LAT, lo in LON
        @test get_ground("epsilonmap", la, lo) == get_epsilon(la, lo)
        @test get_ground("sigmamap", la, lo) == get_sigma(la, lo)
        @test get_ground("codemap", la, lo) == get_groundcode(la, lo)

        @test get_epsilon(la, lo) ≈ GROUND[get_groundcode(la, lo)].ϵᵣ
        @test get_sigma(la, lo) ≈ GROUND[get_groundcode(la, lo)].σ atol=1e-8
    end

    tx = TRANSMITTER[:NAA]
    rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())
    grounds, distances = groundsegments(tx.latitude, tx.longitude, rx.latitude, rx.longitude)
    @test minimum(diff(distances)) >= 20e3
    @test length(grounds) != length(unique(grounds))
    @test length(grounds) == length(distances)

    grounds, distances = groundsegments(tx.latitude, tx.longitude, rx.latitude, rx.longitude;
                                        require_permittivity_change=true)
    @test length(grounds) == length(unique(grounds))
    @test length(grounds) == length(distances)

    # Magnetic field

    az, _, dist, _ = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude)

    # Based on https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
    declination = -16.0378
    calcb = BField(50.6701e-6, deg2rad(67.7737), deg2rad(az-declination))

    bfield = igrf(az, tx.latitude, tx.longitude, 2020)
    all(isapprox(getfield(bfield,f), getfield(calcb,f); rtol=1e-4) for f in fieldnames(BField))
end

###
# Plots
#==

# Ground

tx = TRANSMITTER[:NAA]
rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())

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
grounds, distances = groundsegments(tx.latitude, tx.longitude, rx.latitude, rx.longitude;
    require_permittivity_change=false)
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

==#