#!/bin/bash

configs=
exp=

until [ -z $1 ]
do
    configs="$configs ${exp%/}"
    exp=$1
    shift
done

exp=$(basename $exp)
exp=${exp%.out}.out

header="             :"
files=
for config in $configs
do
    wtrail="${config}                             "
    header="${header} ${wtrail:0:30}:"
    files="$files ${config}/$exp"
done

for file in $files
do
    cat $file | grep TIME | colrm 1 5 | awk '{printf("%-45s\n", $0) }' > $file.times
done

A=
for file in $files
do
    if [ -z $A ]
    then
        A=$file.times
    else
        join -t: $A $file.times > tmp.times
        mv tmp.times joined.times
        A=joined.times
    fi
done

mv joined.times tmp.times
echo "$header" > joined.times
cat tmp.times >> joined.times
rm tmp.times

echo '================================================================'
echo 'cat joined.times'
echo '================================================================'
cat joined.times
