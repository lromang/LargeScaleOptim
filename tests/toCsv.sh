#! /bin/bash

file=$(echo $1 | awk -F '.' '{print $1}')

cat $1 | grep -oE explored.* | awk -F ';' '{print $1 , $2}' | sed -e 's/[a-z]*//g' | sed -r 's/ +//g' | sed -r 's/^=//g' | sed 's/=/,/g' | grep -v '^$'> $file.csv

