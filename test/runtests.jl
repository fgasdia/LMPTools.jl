using Test, LMPTools
using Dates
using GeographicLib
using LongwaveModePropagator
const LMP = LongwaveModePropagator

const GROUND_DATA = LMPTools.GROUND_DATA
const GROUND_LAT = GROUND_DATA["lat"]
const GROUND_LON = GROUND_DATA["lon"]

function bfieldlatlons(tx, rx)
    dists = 0:100e3:range(tx, rx)
    line = GeodesicLine(tx, rx)
    
    lats = Vector{Float64}(undef, length(dists))
    lons = similar(lats)
    for i in eachindex(dists)
        for (i, d) in enumerate(dists)
            lo, la, _ = forward(line, d)
            lats[i] = la
            lons[i] = lo
        end
    end
    return dists, lats, lons
end

@testset "LMPTools" begin
    tx = TRANSMITTER[:NAA]
    rx = Receiver("Boulder", 40.01, -105.244, 0.0, VerticalDipole())

    # Extended functions
    gl1 = GeodesicLine(tx, rx)
    gl2 = GeodesicLine(tx.longitude, tx.latitude; lon2=rx.longitude, lat2=rx.latitude)
    for f in fieldnames(GeodesicLine)
        # testing with ≈ because some values are ~1e-315
        @test getfield(gl1, f) ≈ getfield(gl2, f)
    end

    # https://www.cqsrg.org/tools/GCDistance/
    @test range(tx, rx) ≈ 3142.0e3 atol=100

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

    for la in GROUND_LAT, lo in GROUND_LON
        @test get_ground("epsilonmap", la, lo) == get_epsilon(la, lo)
        @test get_ground("sigmamap", la, lo) == get_sigma(la, lo)
        @test get_ground("codemap", la, lo) == get_groundcode(la, lo)

        @test get_epsilon(la, lo) ≈ GROUND[get_groundcode(la, lo)].ϵᵣ
        @test get_sigma(la, lo) ≈ GROUND[get_groundcode(la, lo)].σ atol=1e-8
    end

    grounds, distances = groundsegments(tx, rx)
    @test minimum(diff(distances)) >= 20e3
    @test length(grounds) != length(unique(grounds))
    @test length(grounds) == length(distances)

    grounds, distances = groundsegments(tx, rx; require_permittivity_change=true)
    @test length(grounds) == length(unique(grounds))
    @test length(grounds) == length(distances)

    # Magnetic field

    az, _, dist, _ = inverse(tx.longitude, tx.latitude, rx.longitude, rx.latitude)

    # Based on https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
    declination = -16.0378
    calcb = BField(50.6701e-6, deg2rad(67.7737), deg2rad(az-declination))

    bfield = igrf(az, tx.latitude, tx.longitude, 2020)
    @test all(isapprox(getfield(bfield,f), getfield(calcb,f); rtol=1e-4) for f in fieldnames(BField))

    bfield2 = only(igrf(tx, rx, 2020, 0.0))
    @test all(isapprox(getfield(bfield2,f), getfield(bfield,f)) for f in fieldnames(BField))

    # chaos does not match igrf, but they should be reasonably close
    bfield3 = chaos(az, tx.latitude, tx.longitude, 2020)
    @test all(isapprox(getfield(bfield3,f), getfield(bfield,f); rtol=0.1) for f in fieldnames(BField))

    ranges, lats, lons = bfieldlatlons(tx, rx)
    chaos_bfields = chaos(tx, rx, 2020, ranges)
    igrf_bfields = igrf(tx, rx, 2020, ranges)

    @test chaos_bfields == chaos(az, lats, lons, 2020)
    for i in 1:length(ranges)
        @test chaos_bfields[i] == chaos(az, lats[i], lons[i], 2020)
        @test all(isapprox(getfield(chaos_bfields[i],f), getfield(igrf_bfields[i],f); rtol=1) for f in fieldnames(BField))
    end

    # Try loading a different CHAOS mat file
    load_CHAOS_matfile(joinpath(LMPTools.project_path("data"), "CHAOS-7.8.mat"))
    bfield4 = chaos(az, tx.latitude, tx.longitude, 2020)
    @test all(isapprox(getfield(bfield4,f), getfield(bfield,f); rtol=0.1) for f in fieldnames(BField))

    load_CHAOS_matfile(joinpath(LMPTools.project_path("data"), "CHAOS-7.11.mat"))
    bfield6 = chaos(az, tx.latitude, tx.longitude, 2020)
    @test all(!isequal(getfield(bfield4,f), getfield(bfield6,f)) for f in fieldnames(BField))

    load_CHAOS_matfile(joinpath(LMPTools.project_path("data"), "CHAOS-8.3.mat"))
    bfield7 = chaos(az, tx.latitude, tx.longitude, 2020)
    @test all(!isequal(getfield(bfield4,f), getfield(bfield7,f)) for f in fieldnames(BField))

    load_CHAOS_matfile(v"7.8")
    bfield5 = chaos(az, tx.latitude, tx.longitude, 2020)
    @test all(isequal(getfield(bfield4,f), getfield(bfield5,f)) for f in fieldnames(BField))

    @test LMPTools.newestchaos() == "CHAOS-8.3.mat"

    # Zenith angle

    dt = DateTime(2021, 2, 1, 12, 45)
    lon, lat = forward(gl1, 1000e3)
    @test zenithangle(lat, lon, dt) == zenithangle(lat, lon, 2021, 2, 1, 12.75)

    # Check out https://ssd.jpl.nasa.gov/horizons.cgi for comparison numbers
    @test zenithangle(lat, lon, dt) ≈ 90 - 0.3106 atol=0.01
    @test zenithangle(lat, lon, 2021, 2, 1, 16) ≈ 90 - 24.9418 atol=0.01
    @test zenithangle(lat, lon, 2021, 2, 1, 21.25) ≈ 90 - 10.5537 atol=0.01

    lon, lat = -110, 62
    @test zenithangle(lat, lon, 2021, 2, 1, 9.5) ≈ 90 - -40.5603 atol=0.01
    @test zenithangle(lat, lon, 2021, 2, 1, 20) ≈ 90 - 10.9651 atol=0.01
    @test zenithangle(lat, lon, 2021, 9, 15, 5) ≈ 90 - -20.1376 atol=0.01
    @test zenithangle(lat, lon, 2021, 9, 15, 19.25) ≈ 90 - 30.7122 atol=0.01

    @test isday(zenithangle(lat, lon, 2021, 2, 1, 16)) == true  # day
    @test isday(zenithangle(lat, lon, 2021, 2, 1, 22)) == true
    @test isday(zenithangle(lat, lon, 2021, 2, 2, 1)) == false  # night
    @test isday(zenithangle(lat, lon, 2021, 2, 2, 6)) == false
    @test isday(zenithangle(lat, lon, 2021, 2, 2, 8)) == false
    @test isday(zenithangle(lat, lon, 2021, 2, 2, 11)) == false

    # Ionospheres

    dt = DateTime(2020, 3, 1, 2)
    lats = 15:89
    lons = -160:-60
    szas = [zenithangle(la, lo, dt) for la in lats, lo in lons]

    # Compare number and array-specialized form
    x = [flatlinearterminator(s) for s in szas]
    hprimes, betas = getindex.(x,1), getindex.(x,2)
    hprimes_arr, betas_arr = flatlinearterminator(szas)
    @test hprimes == hprimes_arr
    @test betas == betas_arr

    h, b = smoothterminator(50)
    @test h ≈ 74  # hp_day
    @test b ≈ 0.3  # b_day

    h, b = smoothterminator(130)
    @test h ≈ 86  # hp_night
    @test b ≈ 0.5  # b_night

    h, b = smoothterminator(95)
    @test 74 < h < 86
    @test 0.3 < b < 0.5

    # Sanity tests on fourierperturbation
    coeff = zeros(5)
    pert = fourierperturbation.(szas, (coeff,))
    @test iszero(pert)

    coeff = [1, 0, 0, 0, 0]
    pert = fourierperturbation.(szas, (coeff,))
    @test all(isequal(1), pert)
end
