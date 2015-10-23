set term gif animate
set output "animate.gif"
stats 'heat_diffusion.dat' nooutput
minT=0
maxT=1.591549
set cbrange [minT:maxT]
do for [i=1:int(STATS_blocks)] {
  plot 'heat_diffusion.dat' index (i-1)  using 1:2:3 with image
}
