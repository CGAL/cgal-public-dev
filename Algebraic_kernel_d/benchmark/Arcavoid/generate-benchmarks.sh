#!/bin/bash

GENERATOR=./build/release/generators

if ( false ); then
    mkdir -p data/rand_u_d
    for n in 0{1..9}0 1{0..9}0 200 500; do
	for L in 0010 0050 0100 0200 0500 1000; do
            for seed in {0..9}; do
		$GENERATOR u mpsolve $n $L $seed > data/rand_u_d/rud-$n-$L-$seed.in
            done
	done
    done
fi

if ( false ); then
    mkdir -p data/rand_u_s
    for n in 0{1..9}0 1{0..9}0 200 500; do
	for L in 0010 0050 0100 0200 0500 1000; do
            for seed in {0..9}; do
		$GENERATOR us mpsolve $n $L $seed > data/rand_u_s/rus-$n-$L-$seed.in
            done
	done
    done
fi

if ( false ); then
    mkdir -p data/e
    for n in 0{2..9}{0,5} 1{0..9}{0,5} 200; do
	$GENERATOR e mpsolve $n 0 0 > data/e/e-$n.in
    done
fi