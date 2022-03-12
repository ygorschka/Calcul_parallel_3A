set autoscale

set xlabel 'cc'

splot 'sol000.txt' w p title 'proc 0' \
, 'sol001.txt' w p title 'proc 1' \
, 'sol002.txt' w p title 'proc 2' \
, 'sol003.txt' w p title 'proc 3'

#set term png
#set output "couette.png"
#replot
#set term x11
