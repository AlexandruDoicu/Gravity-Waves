set pm3d map interpolate 0,0
set palette gray positive gamma 1.25
set term png
set output "WaveResultsWP/SrcFctWP.png"
splot "WaveResultsWP/SrcFctWP.dat" using 1:2:3
