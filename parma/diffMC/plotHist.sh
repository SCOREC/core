#!/bin/bash
gnuplot << EOF
set term png size 1600, 1200
set output "$1.png"
set style data histograms
set style histogram cluster
set style fill solid 1.0 border lt -1
set yrange [0:*]
p '$1' u 2:xticlabels(1) t 'bin weight'
EOF
