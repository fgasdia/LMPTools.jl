using Dates, Printf
using Plots
using CSV, DataFrames
using LongwaveModePropagator
using LMPTools

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

# Zenith angle

dt = DateTime(2021, 2, 1)
lats = 15:89
lons = -160:-60
szas = [zenithangle(la, lo, dt) for la in lats, lo in lons]
heatmap(lons, lats, szas,
        color=cgrad(:starrynight, [50, 70, 80, 85, 90, 95, 100, 110, 130]/180, rev=true), clims=(0, 180))



function gpmap(file, lats::AbstractRange, lons::AbstractRange, data)
    size(data) == (length(lats), length(lons)) || throw(ArgumentError("data not compatible with lats, lons"))
    issorted(lons) && issorted(lats) || throw(ArgumentError("lats and lons must be sorted"))

    δlat = step(lats)/2
    δlon = step(lons)/2

    open(file, "w") do f
        for j in eachindex(lons)
            for i in eachindex(lats)
                
                str = @sprintf("%f, %f, %f\n", lons[j]-δlon, lats[i]-δlat, szas[i,j])
                write(f, str)
                str = @sprintf("%f, %f, %f\n", lons[j]-δlon, lats[i]+δlat, szas[i,j])
                write(f, str)
                str = @sprintf("%f, %f, %f\n", lons[j]+δlon, lats[i]+δlat, szas[i,j])
                write(f, str)
                str = @sprintf("%f, %f, %f\n", lons[j]+δlon, lats[i]-δlat, szas[i,j])
                write(f, str)
                
                write(f, "\n")
            end
        end
    end
end

dt = DateTime(2021, 2, 1)
lats = 15:89
lons = -150:-65
szas = [zenithangle(la, lo, dt) for la in lats, lo in lons]
gpmap("szas.csv", lats, lons, szas)

