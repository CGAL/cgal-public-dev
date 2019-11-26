#!/bin/bash

if [ $# -ne 1 ]
then
    echo "usage: generate-times fic"
    echo " Which extracts results for homotopy tests from file 'fic'"
    echo " and outputs the data in the two files computation-time-reduce-surface.dat"
    echo " and computation-time-path-homotopy.dat."
    exit 1
fi

ficIn="$1"

# 1) Extract computation times of homotopy tests.

# 1.1) Extract the different info.
echo "Path1" > restmp1.txt
grep "Random seed: " "${ficIn}" | cut -d ' ' -f 6 >> restmp1.txt
echo "Path2" > restmp2.txt
grep "Random seed: " "${ficIn}" | cut -d ' ' -f 12 >> restmp2.txt
echo "TimeContractible" > restmp3.txt
grep "\[TIME\] is_contractible: " "${ficIn}" | cut -d ' ' -f 3 >> restmp3.txt
echo "TimeHomotopy" > restmp4.txt
grep "\[TIME\] are_freely_homotopic: " "${ficIn}" | cut -d ' ' -f 3 >> restmp4.txt

# 1.2) Regroup the different info in different columns of a same file.
echo "# Size_of_path1; Size_of_path 2; time_of_contractible_test; time_of_homotopy_test." > computation-time-path-homotopy.dat
paste -d '\t' restmp1.txt restmp2.txt restmp3.txt restmp4.txt >> computation-time-path-homotopy.dat

# 2) Extract 2-map size and simplification times

echo "# ==========#INITIAL-MAP==========     ==========#REDUCED-MAP========" > computation-time-reduce-surface.dat
echo "#darts   #vertices #edges   #faces     #darts #vertices #edges #faces    GlobalTime" >> computation-time-reduce-surface.dat

# 2.1) Extract the different info.
grep "Initial map:" "${ficIn}" | cut -d '=' -f 2 | cut -d ',' -f 1 > restmp1.txt
grep "Initial map:" "${ficIn}" | cut -d '=' -f 3 | cut -d ',' -f 1 > restmp2.txt
grep "Initial map:" "${ficIn}" | cut -d '=' -f 4 | cut -d ',' -f 1 > restmp3.txt
grep "Initial map:" "${ficIn}" | cut -d '=' -f 5 | cut -d ',' -f 1 > restmp4.txt

grep "Reduced map:" "${ficIn}" | cut -d '=' -f 2 | cut -d ',' -f 1 > restmp5.txt
grep "Reduced map:" "${ficIn}" | cut -d '=' -f 3 | cut -d ',' -f 1 > restmp6.txt
grep "Reduced map:" "${ficIn}" | cut -d '=' -f 4 | cut -d ',' -f 1 > restmp7.txt
grep "Reduced map:" "${ficIn}" | cut -d '=' -f 5 | cut -d ',' -f 1 > restmp8.txt

grep "\[TIME\] Total time for computation of reduced map: " "${ficIn}" | cut -d ' ' -f 9  > restmp9.txt

# 2.2) Regroup the different info in different columns of a same file.
paste -d '\t' restmp1.txt restmp2.txt restmp3.txt restmp4.txt restmp5.txt restmp6.txt restmp7.txt restmp8.txt restmp9.txt >> computation-time-reduce-surface.dat

# 3) Erase temp files.
rm -f restmp1.txt restmp2.txt restmp3.txt restmp4.txt restmp5.txt restmp6.txt restmp7.txt restmp8.txt restmp9.txt
