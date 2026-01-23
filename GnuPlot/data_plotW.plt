set title "Horizontal Velocity [m/s]"
set pm3d map interpolate 0,0
set xlabel "Time [hr]"
set ylabel "Altitude [km]" offset-2.0
#set xrange [0:24]
#set yrange [0:24]
#set zrange [0:100]
#set cbrange [-100:100]
set lmargin 10
set palette gray positive gamma 1.25
set term png
set output "WaveResultsWP/WVelocityWP.png"
splot "WaveResultsWP/WVelocityWP.dat" using 1:2:3
