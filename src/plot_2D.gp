unset key
set palette defined ( 0 "blue", 0.5 "white", 1 "red" )
plot "sol000.txt" u 1:2:3 w p pt 7 ps 3 lc palette z \
, "sol001.txt" u 1:2:3 w p pt 7 ps 3 lc palette z \
, "sol002.txt" u 1:2:3 w p pt 7 ps 3 lc palette z \
, "sol003.txt" u 1:2:3 w p pt 7 ps 3 lc palette z
