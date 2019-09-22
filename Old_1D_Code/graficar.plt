load 'gnuplot_config.ini'
set title "Concentraci{\363}n a lo largo del reactor para diferentes tiempos \n usando un esquema implicito" font ",25"
p "results.dat" i 0:0 u 1:2 w lp ls 1 t "t=0.0s",\
"" i 100:100 w lp ls 2 t "t=0.1s",\
"" i 200:200 w lp ls 3 t "t=0.2s",\
"" i 400:400 w lp ls 4 t "t=0.4s",\
"" i 800:800 w lp ls 5 t "t=0.8s",\
"" i 1600:1600 w lp ls 6 t "t=1.6s",\
"" i 3200:3200 w l ls 7 t "t=3.2s"
