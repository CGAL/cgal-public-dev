#!/bin/bash
for filename in data/*; do
    ./poisson_smooth_test "$filename"
    echo -e "\n \n"
done
