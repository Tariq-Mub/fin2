set term png
set o "mc-pi1.png"

set logscale x 10
set logscale y 10

p "mc_pi.res" u 1:3,  1/sqrt(x)
