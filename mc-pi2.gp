unset key

set term png
set o "mc-pi2.png"

set logscale x 10

p [700000:800000000] [3.139:3.144] "mc_pi.res" w erro, pi
