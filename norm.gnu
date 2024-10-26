#!/usr/bin/gnuplot -persist

set term pdf
set output '3D.pdf'
set hidden3d
set pm3d
unset surface
#set view 50,220
set xrange [__FILL_:__FILL_]
set yrange [__FILL_:__FILL_]
unset ztics
set view 0,0
set contour base
set cntrparam order 4
set cntrparam linear
set cntrparam levels discrete 0.7
set cntrparam points 5
splot 'fort.8' u 1:2:3 notitle w l


