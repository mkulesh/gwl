set pm3d map
set palette color negative
set xtics 0,8,128
unset ytics 

set multiplot

set size 0.856,0.32
set origin 0.052,0
plot 'exmDeltaSigASC.dat' w l notitle

set size 1.1,0.55
set origin -0.07,0.18
splot 'exmDeltaCwtGnu.dat' u 1:2:4 w pm3d notitle

set size 1.1,0.55
set origin -0.07,0.5
splot 'exmDeltaCwtGnu.dat' u 1:2:3 w pm3d notitle

unset multiplot


