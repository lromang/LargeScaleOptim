#! /bin/bash

dir=$(pwd)

for i in $(ls | grep .*.txt); do ../toCsv.sh "$dir/$i"; done

