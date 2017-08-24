#!/bin/bash

for file in tests/*
do
  echo "$file"
  ./main "$file"
done 
