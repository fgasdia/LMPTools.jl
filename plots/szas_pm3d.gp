set terminal pngcairo size 600,450 lw 1 fontscale 1 #truecolor transparent
set output 'szas_pm3d.png'

set loadpath 'C:\\Users\\forrest\\gnuplot\\' 'C:\\Users\\forrest\\gnuplot\\maps\\'

load 'common.cfg'
load 'mydistinctcolors.pal'

set angle degrees

unset mouse  # tracking doesn't report longitude, latitude

d2r(x) = x*pi/180.0
cot(x) = 1.0/tan(x)
Re = 1.0

ϕ1 = 40.0  # standard parallel
λ0 = -102.0  # reference longitude (centered at 0, 0)
ϕ0 = 0.0  # reference latitude

# Lambert conformal conic projection (single parallel form, in degrees!)
n = sin(ϕ1)
F = cos(ϕ1)*tan(45.0 + 0.5*ϕ1)**n/n
θ(λ) = n*(λ - λ0)
ρ(ϕ) = Re*F*cot(45.0 + 0.5*ϕ)**n
ρ0 = ρ(ϕ0)

x_lambert(λ, ϕ) = ρ(ϕ)*sin(θ(λ))
y_lambert(λ, ϕ) = ρ0 - ρ(ϕ)*cos(θ(λ))

unset xtics
unset ytics
unset border

set bmargin at screen 0.05
set tmargin at screen 0.95
set lmargin at screen 0.01
set rmargin at screen 0.85

ymin = d2r(39)
ymax = d2r(78)
xhalfwidth = d2r(35)  # +/- this distance from λ0

set size ratio (ymax - ymin)/(2*xhalfwidth) 1,1
# set size ratio 0.5 1,1

# ranges are sort of meaningless because the map is now curved, but scaled approximately to radians
set xrange [-xhalfwidth:xhalfwidth]
set yrange [ymin:ymax]
set zrange [0:1]
set cbrange [0:180]

set colorbox user size 0.03,0.5 orig 0.87,0.25

# set style fill transparent solid 0.8 noborder

set view map
set hidden3d front  # for some reason, this command results in reverse z ordering, so higher z is plotted below lower z

splot 'szas_3col.csv' u (x_lambert($1,$2)):(y_lambert($1,$2)):(0):3 w pm3d notitle, \
	  for [λ=-180:180:10] [ϕ=30:90] '+' u (x_lambert(λ,ϕ)):(y_lambert(λ,ϕ)):(1) w l lc rgb '#dddddd' lw 0.5 notitle, \
	  for [ϕ=30:90:10] [λ=-180:180:2] '+' u (x_lambert(λ,ϕ)):(y_lambert(λ,ϕ)):(1) w l lc rgb '#dddddd' lw 0.5 notitle, \
	  'ne_50m_coastline_gp.csv' u (x_lambert($1,$2)):(y_lambert($1,$2)):(0) w l lc 'black' lw 0.5 notitle, \
	  'ne_50m_admin_1_states_provinces_lakes_gp.csv' u (x_lambert($1,$2)):(y_lambert($1,$2)):(0) w l lc 'black' lw 0.5 notitle, \
	 