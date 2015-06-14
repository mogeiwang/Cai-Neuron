set term pdfcairo font "Times New Roman, 8"
set output "plot1-4.pdf"
set multiplot layout 4,1
plot "vol0.txt" w l
plot "vol1.txt" w l
plot "vol2.txt" w l
plot "vol3.txt" w l
set label 1 "time (S)" at 3.0,-0.1
set label 2 "voltage of neuron #i" at -0.3,1.2
unset multiplot
