set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'convergencia_erro_log.png'

# Label the axes as log values
set xlabel 'log(h)'
set ylabel 'log(erro absoluto)'

# Define a linear function for the transformed data
f(x) = inclinacao*x + berro

# Fit the function using log-transformed data
fit f(x) '/home/sagan/SL_NEWTON_RAPHSON/ORDEM/convergencia_erro.dat' using (log($1)):(log($2)) via inclinacao,berro

# Print the values of a and b on the graph
set label sprintf("a = %g", inclinacao) at graph 0.1, 0.9
set label sprintf("b = %g", berro) at graph 0.1, 0.85

# Plot the transformed data and the fitted curve
plot '/home/sagan/SL_NEWTON_RAPHSON/ORDEM/convergencia_erro.dat' using (log($1)):(log($2)) with points pt 7 ps 1.5 title 'Log Data', \
    f(x) with lines lw 2 title 'Fit escala logartimica'