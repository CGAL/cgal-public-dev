#!/bin/bash
  for filename in data/*; do
    echo -e "$filename"
    for i in $(seq 0.0001 0.00099 0.01); do
      ./poisson_smooth_test "$filename" -approx "$i" -sm_radius 0.5
      echo -e "\n"
    done
    echo -e "*****************"
  done
  echo -e "******************************\n \n******************************\n \n"

  for i in $(seq -0.5 0.1 0.5); do
    ./poisson_smooth_test data/sphere_hole.pwn -isovalue "$i"
    echo -e "\n \n"
  done
