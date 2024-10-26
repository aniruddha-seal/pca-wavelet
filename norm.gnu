#!/usr/bin/gnuplot -persist

set term pdf
set output '3D.pdf'
set hidden3d
set pm3d
unset surface
#set view 50,220
set xrange [0:4.0]
set yrange [700:4000]
unset ztics
set view 0,0
set contour base
set cntrparam order 4
set cntrparam linear
set cntrparam levels discrete 0.7
set cntrparam points 5
splot 'fort.8' u 1:2:3 notitle w l


