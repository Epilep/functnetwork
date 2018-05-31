set terminal pngcairo  transparent enhanced font "arial,12" #fontscale 1.0 size 600, 400 
set output 'functional_matrix_old.png'
unset key
#set view map scale 1
set style data lines
set xtics border in scale 0,0 mirror norotate  #autojustify
set ytics border in scale 0,0 mirror norotate  #autojustify
set ztics border in scale 0,0 nomirror norotate  #autojustify
#unset cbtics
#set rtics axis in scale 0,0 nomirror norotate  autojustify
#set title "Heat Map generated by 'plot' from a stream of XYZ values\nNB: Rows must be separated by blank lines!" 
set xrange [ 0 : 67 ] noreverse nowriteback
set yrange [ 0 : 67 ] noreverse nowriteback
set cblabel "Association" 
set cbrange [ 0.00000 : 0.350000 ] noreverse nowriteback
#set palette rgbformulae -7, 2, -7
set palette rgb 34,35,36
#DEBUG_TERM_HTIC = 119
#DEBUG_TERM_VTIC = 119
## Last datafile plotted: "$map2"
plot 'functional_matrix_old.dat' using 2:1:3 with image
