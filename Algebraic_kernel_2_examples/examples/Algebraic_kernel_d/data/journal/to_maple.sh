#!/bin/bash

out="maple_in"

echo 'F :=' > $out
echo '[' >> $out

for i in $*; do
    echo '  [' >> $out
    echo '    [' >> $out
    cat $i | sed -E -e 's/^(.*)$/      &,/' >> $out
    echo '      0' >> $out
    echo '    ],' >> $out
    echo '  "'$i'"' >> $out
    echo '  ],' >> $out
done

echo '  0' >> $out
echo ']:' >> $out
