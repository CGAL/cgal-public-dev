#!/bin/bash

TIME="gtime -v"

BUILD_DIR=$HOME/CGAL/branches/Arcavoid/Arcavoid-akobel/examples/build/release

CONVERT_POLYNOMIAL=$BUILD_DIR/convert_polynomial

ABSOLUTE_VOID=$BUILD_DIR/Absolute_void
BITSTREAM_DESCARTES=$BUILD_DIR/Bitstream_descartes
DESCARTES=$BUILD_DIR/Descartes
RS=$BUILD_DIR/RS
SAGE=$HOME/bin/sage
UNISOLVE="$HOME/archives/MPSolve-2.2/unisolve -o10000 -d1 -Ob -Db"

GNUPLOT=/usr/bin/gnuplot
GNUPLOT=/usr/local/bin/gnuplot

loops=1
instances=""
complex=false
cleanup=false
latex=false
TIMEOUT=

SOLVER_A=true   # Absolute_void
SOLVER_B=true   # Bitstream_descartes
SOLVER_D=true   # Descartes
SOLVER_M=false  # Maple
SOLVER_P=false  # Pari/GP
SOLVER_R=true   # RS
SOLVER_S=true   # Sage
SOLVER_U=true   # MPSolve (unisolve)

while getopts "n:lLkKt:abdmprusABDMPRUS" flag
do
    if [ $flag == "n" ]; then loops=$OPTARG; fi

    if [ $flag == "l" ]; then latex=true; fi
    if [ $flag == "L" ]; then latex=false; fi

    if [ $flag == "k" ]; then cleanup=false; fi
    if [ $flag == "K" ]; then cleanup=true; fi

    if [ $flag == "c" ]; then complex=true; fi
    if [ $flag == "C" ]; then complex=false; fi

    if [ $flag == "t" ]; then
	MAX_TIME="$OPTARG"
	TIMEOUT="timeout $MAX_TIME"
    fi

    if [ $flag == "a" ]; then SOLVER_A=true; fi
    if [ $flag == "b" ]; then SOLVER_B=true; fi
    if [ $flag == "d" ]; then SOLVER_D=true; fi
    if [ $flag == "m" ]; then SOLVER_M=true; fi
    if [ $flag == "p" ]; then SOLVER_P=true; fi
    if [ $flag == "r" ]; then SOLVER_R=true; fi
    if [ $flag == "s" ]; then SOLVER_S=true; fi
    if [ $flag == "u" ]; then SOLVER_U=true; fi

    if [ $flag == "A" ]; then SOLVER_A=false; fi
    if [ $flag == "B" ]; then SOLVER_B=false; fi
    if [ $flag == "D" ]; then SOLVER_D=false; fi
    if [ $flag == "M" ]; then SOLVER_M=false; fi
    if [ $flag == "P" ]; then SOLVER_P=false; fi
    if [ $flag == "R" ]; then SOLVER_R=false; fi
    if [ $flag == "S" ]; then SOLVER_S=false; fi
    if [ $flag == "U" ]; then SOLVER_U=false; fi
done

echo -n "Benchmarking:"
if [ $SOLVER_A == true ]; then echo -n " Absolute_void"; fi
if [ $SOLVER_B == true ]; then echo -n " Bitstream_descartes"; fi
if [ $SOLVER_D == true ]; then echo -n " Descartes"; fi
if [ $SOLVER_M == true ]; then echo -n " Maple"; fi
if [ $SOLVER_P == true ]; then echo -n " Pari/GP"; fi
if [ $SOLVER_R == true ]; then echo -n " RS"; fi
if [ $SOLVER_S == true ]; then echo -n " Sage"; fi
if [ $SOLVER_U == true ]; then echo -n " MPSolve"; fi
echo

shift $((OPTIND-1))

NOW=`date +%Y_%m_%d__%H_%M`
PWD=`pwd`
ID=
if [ -d "$1" ]
then
    ID="$(echo "${1%/}" | sed -e 's|/|_|g' -e 's| |_|g')"
else
    ID="$(echo "$(basename "${1%/}")" | sed -e 's|/|_|g' -e 's| |_|g')"
    ID="${ID%.in}"
fi
ID="${ID#data_}"
WDIR="$PWD/exp_results__${ID}__${NOW}__$(printf "%05i" $$)"
mkdir $WDIR
RESULTS=$WDIR/results.tbl

echo "Output directory: $WDIR"
if [ $cleanup == true ]; then
    echo "(will be removed after benchmarks are run)"
fi

echo

SEP=""
EOL=""
if [ $latex == true ]; then
    SEP="&"
    EOL='\\'
fi
DATAFILE_LENGTH=40
DEG_BITS_LENGTH=5
INSTANCE_LENGTH=10

format_datafile () {
    echo -n "$(printf "%-${DATAFILE_LENGTH}s" "$1")" | tee -a $RESULTS
}
format_degree_bitsize () {
    echo -n " ${SEP} $(printf "%${DEG_BITS_LENGTH}s" $1)" | tee -a $RESULTS
    echo -n " ${SEP} $(printf "%${DEG_BITS_LENGTH}s" $2)" | tee -a $RESULTS
}
format_instance () {
#    entry="          $(printf "%.02f" ${1})"
#    echo -n " ${SEP}${entry: -10}" | tee -a $RESULTS
    if [ "$1" = ">$MAX_TIME" ]; then
	echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" "$1")" | tee -a $RESULTS
    else
	echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" "$(printf "%.02f" "$1")")" | tee -a $RESULTS
    fi
}
finish_line () {
    echo "${EOL}" | tee -a $RESULTS
}

convert_polynomial () {
    ./build/release/convert_polynomial $1
}

run_Absolute_void () {
    mkdir -p absolute_void
    $TIMEOUT $TIME -o absolute_void/$1.time \
        $ABSOLUTE_VOID _input/$1.mpsolve &> absolute_void/$1.out
#    format_instance $(cat absolute_void/$1.time | grep User | cut -d: -f2)
    if [ $? -gt 123 ]; then
	format_instance ">$MAX_TIME"
	return
    fi
    format_instance $(cat absolute_void/$1.out | grep "TIME FOR ISOLATION" | cut -d: -f2)
    
    cat absolute_void/$1.out | grep -e "complex root [0-9]*:" | cut -d"(" -f3 | cut -d";" -f1 | sed -e "s/\+ I\*/\t/" > absolute_void/$1.gplot
}

run_Bitstream_descartes () {
    mkdir -p bitstream_descartes
    $TIMEOUT $TIME -o bitstream_descartes/$1.time \
        $BITSTREAM_DESCARTES _input/$1.cgal &> bitstream_descartes/$1.out
    if [ $? -gt 123 ]; then
	format_instance ">$MAX_TIME"
	return
    fi
    format_instance $(cat bitstream_descartes/$1.out | grep Time | cut -d: -f2 | cut -d" " -f2)
}

run_Descartes () {
    mkdir -p descartes
    $TIMEOUT $TIME -o descartes/$1.time \
        $DESCARTES _input/$1.mpsolve &> descartes/$1.out
    if [ $? -gt 123 ]; then
	      format_instance ">$MAX_TIME"
	      return
    fi
    format_instance $(cat descartes/$1.out | grep "TIME FOR ISOLATION" | cut -d: -f2)
}

run_RS () {
    mkdir -p rs
    $TIMEOUT $TIME -o rs/$1.time \
        $RS _input/$1.cgal &> rs/$1.out
    if [ $? -gt 123 ]; then
	      format_instance ">$MAX_TIME"
	      return
    fi
    format_instance $(cat rs/$1.out | grep "TIME FOR ISOLATION" | cut -d: -f2)
}

run_maple () {
    format_instance "-1"
}

run_pari () {
    format_instance "-1"
}

#if [ $SOLVER_S == true ]; then
#    coproc sage_server { $SAGE; }
#    echo 'R.<x> = QQ[];' >&${sage_server[1]}
#    echo 'from sage.rings.polynomial.real_roots import *;' >&${sage_server[1]}
#    echo 'from sage.rings.polynomial.complex_roots import *;' >&${sage_server[1]}
#fi

run_sage_real () {
    mkdir -p sage_real
    $TIMEOUT $TIME -o sage_real/$1.time \
        $SAGE < _input/$1.sage_real &> sage_real/$1.out
    if [ $? -gt 123 ]; then
	format_instance ">$MAX_TIME"
	return
    fi

#    echo "$(cat _input/$1.sage_real | grep 'f = ')" >&${sage_server[1]}
    
#    echo '%time sols = real_roots (f)' >&${sage_server[1]}
#    echo 'for sol in sols :' >&${sage_server[1]}
#    echo -e '\t["REAL_ROOT", sol[0]]' >&${sage_server[1]}
#    echo -e '' >&${sage_server[1]}
#    echo '"FINISHED"' >&${sage_server[1]}
    
#    REPLY="notfinish"
#    until ( echo "$REPLY" | grep "FINISHED" &> /dev/null ); do
#        read -u ${sage_server[0]}
#        echo "$REPLY" >> sage_real/$1.out
#    done

    format_instance $(cat sage_real/$1.out | grep "CPU times" | sed -e 's/sage://g' | cut -d: -f2 | cut -d" " -f3)
}
run_sage_complex () {
    mkdir -p sage_complex
    $TIMEOUT $TIME -o sage_complex/$1.time \
        $SAGE < _input/$1.sage_complex &> sage_complex/$1.out
    if [ $? -gt 123 ]; then
	format_instance ">$MAX_TIME"
	return
    fi

#    echo "$(cat _input/$1.sage_complex | grep 'f = ')" >&${sage_server[1]}

#    echo '%time sols = complex_roots (f, skip_squarefree=true);' >&${sage_server[1]}
#    echo 'for sol in sols :' >&${sage_server[1]}
#    echo -e '\t["ROOT", real (sol[0]), imag (sol[0])]' >&${sage_server[1]}
#    echo -e '' >&${sage_server[1]}
#    echo '"FINISHED"' >&${sage_server[1]}

#    REPLY="notfinish"
#    until ( echo "$REPLY" | grep "FINISHED" &> /dev/null ); do
#        read -u ${sage_server[0]}
#        echo "$REPLY" >> sage_complex/$1.out
#    done

    format_instance $(cat sage_complex/$1.out | grep "CPU times" | sed -e 's/sage://g' | cut -d: -f2 | cut -d" " -f3)

    cat sage_complex/$1.out | grep "ROOT" | cut -d'[' -f2 | cut -d']' -f1 | sed -e "s/'ROOT', //g" -e 's/\?//g' -e 's/,/ /' > sage_complex/$1.gplot
}

run_unisolve () {
    mkdir -p unisolve
    SEED="$(date +"%N")"
    $TIMEOUT $TIME -o unisolve/$1.time \
        $UNISOLVE "-H${SEED}" _input/$1.mpsolve &> unisolve/$1.out
    if [ $? -gt 123 ]; then
	format_instance ">$MAX_TIME"
	return
    fi
    echo "Seed: $SEED" >> unisolve/$1.out
    format_instance $(cat unisolve/$1.time | grep User | cut -d: -f2)

    cat unisolve/$1.out | grep -e "root [0-9]*" | cut -d"(" -f2 | cut -d")" -f1 | sed -e "s/,/\t/" > unisolve/$1.gplot
}

echo -n "$(printf "%-${DATAFILE_LENGTH}s" Instance)" | tee $RESULTS
echo -n "$(printf " ${SEP} %${DEG_BITS_LENGTH}s" deg)" | tee -a $RESULTS
echo -n "$(printf " ${SEP} %${DEG_BITS_LENGTH}s" bits)" | tee -a $RESULTS
if [ $SOLVER_A == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" Abs_void)" | tee -a $RESULTS; fi
if [ $SOLVER_B == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" Bitstream)" | tee -a $RESULTS; fi
if [ $SOLVER_D == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" Descartes)" | tee -a $RESULTS; fi
if [ $SOLVER_M == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" Maple)" | tee -a $RESULTS; fi
if [ $SOLVER_P == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" Pari/GP)" | tee -a $RESULTS; fi
if [ $SOLVER_R == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" RS)" | tee -a $RESULTS; fi
if [ $SOLVER_S == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" "Sage(RR)")" | tee -a $RESULTS; fi
if [ $SOLVER_S == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" "Sage(CC)")" | tee -a $RESULTS; fi
if [ $SOLVER_U == true ]; then echo -n " ${SEP} $(printf "%${INSTANCE_LENGTH}s" MPSolve)" | tee -a $RESULTS; fi
finish_line

#LINE_LENGTH=$(($(cat $RESULTS | wc -c | cut -d" " -f1) - 1))
#printf "%*s\n" $LINE_LENGTH " " | tr " " "=" | tee -a $RESULTS

mkdir -p ${WDIR}
touch ${WDIR}/INCOMPLETE

for file in $(find $1 -name '*.in' | sort)
do
    sanitized=$(echo $file | sed -e 's|/|_|g' -e 's| |_|g')
    sanitized=${sanitized#data_}
    mkdir -p ${WDIR}/_input
    cp $file ${WDIR}/_input/${sanitized}
    stats="$($CONVERT_POLYNOMIAL ${WDIR}/_input/${sanitized})"
    basename=${sanitized%.in}

    pushd $WDIR &> /dev/null
    
#    ulimit -t 30

    format_datafile $sanitized
    format_degree_bitsize $(echo "$stats" | grep Degree | cut -d: -f2) $(echo "$stats" | grep Bitsize | cut -d: -f2)
    if [ $SOLVER_A == true ]; then run_Absolute_void $basename; fi
    if [ $SOLVER_B == true ]; then run_Bitstream_descartes $basename; fi
    if [ $SOLVER_D == true ]; then run_Descartes $basename; fi
    if [ $SOLVER_M == true ]; then run_maple $basename; fi
    if [ $SOLVER_P == true ]; then run_pari $basename; fi
    if [ $SOLVER_R == true ]; then run_RS $basename; fi
    if [ $SOLVER_S == true ]; then run_sage_real $basename; fi
    if [ $SOLVER_S == true ]; then run_sage_complex $basename; fi
    if [ $SOLVER_U == true ]; then run_unisolve $basename; fi
    finish_line

    mkdir -p gnuplot
    echo "#!$GNUPLOT -p" > gnuplot/$basename.gplot
    echo "plot 0 ls 0;" >> gnuplot/$basename.gplot
    if [ $SOLVER_U == true ]; then echo "replot '$WDIR/unisolve/$basename.gplot' using 1:2 title \"MPSolve\" w p lc 1 pt 4;" >> gnuplot/$basename.gplot; fi
    if [ $SOLVER_S == true ]; then echo "replot '$WDIR/sage_complex/$basename.gplot' using 1:2 title \"Sage\" w p lc 2 pt 2;" >> gnuplot/$basename.gplot; fi
    if [ $SOLVER_A == true ]; then echo "replot '$WDIR/absolute_void/$basename.gplot' using 1:2 title \"Absolute_void\" w p lc 3 pt 1;" >> gnuplot/$basename.gplot; fi
    echo "quit" >> gnuplot/$basename.gplot
    chmod a+x gnuplot/$basename.gplot

    popd &> /dev/null
done

#if [ $SOLVER_S == true ]; then
#    kill ${sage_server_PID}
#fi

mv ${WDIR}/INCOMPLETE ${WDIR}/COMPLETE

if [ $cleanup == true ]; then
    rm -r "$WDIR"
fi
