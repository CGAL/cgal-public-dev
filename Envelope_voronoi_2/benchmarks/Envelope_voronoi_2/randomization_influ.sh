#!/bin/bash
for i in `seq 1 100000`; do
    
    let n=i*100
    
    file_name=dani_points_sorted-$n.in
    echo $file_name
    ./generate/generate -t d -s -p $n > data/points/$file_name
    ./linux_gnu_4.0.4_voronoi_fone/benchVoronoiFone -P data/points -p 0 --type e -s 2 $file_name >> results/randomization_influence_on_sorted_points.txt

done