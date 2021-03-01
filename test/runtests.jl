using Test, LMPTools
using GeographicLib
using LongwaveModePropagator

const GROUND_DATA = LMPTools.GROUND_DATA
const LAT = GROUND_DATA["lat"]
const LON = GROUND_DATA["lon"]


@testset "LMPTools" begin

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
end

###
# Plots
#==

tx = TRANSMITTER[:NAA]
rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())

sigmas = [get_sigma(y, x) for y in LAT, x in LON]
heatmap(LON, LAT, sigmas, clims=(0, 0.11))

epsilons = [get_epsilon(y, x) for y in LAT, x in LON]
heatmap(LON, LAT, epsilons, clims=(0, 20))

pts = [forward(tx.longitude, tx.latitude, az, d) for d in distances]
heatmap(LON, LAT, sigmas, clims=(0, 0.11), xlims=(-110, -65), ylims=(35, 48))
scatter!(getfield.(pts, :lon), getfield.(pts, :lat))

==#