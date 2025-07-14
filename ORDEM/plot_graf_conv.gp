# Set the output to a PNG file (no on-screen display)
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'convergencia_y.png'

# Label the axes
set xlabel 'h'
set ylabel 'erro absoluto'

# Define a 4th order polynomial function
f(x) = a*x**4 + b*x**3 + c*x**2 + d*x + e

# Perform the fit ignoring comment lines (lines starting with #)
fit f(x) '/home/sagan/SL_NEWTON_RAPHSON/ORDEM/convergencia_y.dat' using 1:2 via a,b,c,d,e

# Plot the original data and the fitted curve
set label 1 sprintf("a = %g", a) at graph 0.1, 0.9
set label 2 sprintf("b = %g", b) at graph 0.1, 0.8
set label 3 sprintf("c = %g", c) at graph 0.1, 0.7
set label 4 sprintf("d = %g", d) at graph 0.1, 0.6
set label 5 sprintf("e = %g", e) at graph 0.1, 0.5
plot '/home/sagan/SL_NEWTON_RAPHSON/ORDEM/convergencia_y.dat' using 1:2 with points pt 7 ps 1.5 title 'Dados', \
     f(x) with lines lw 2 title 'Fit de ordem 4'
