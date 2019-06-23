#!/bin/bash

folderpath="/home/fatemeh/Leishmania_2019/Leish_paper_first_round_input_data/"
filenames=$(ls -p "${folderpath}"  | grep -i .eps$) #-v / | grep "*.eps") #| awk -F" " '{print $9}')

#printf "%s\n" "${IDS[@]}" | sort -u
declare -A setA
for name in $filenames; do 
    SUBSTRING=$(echo $name| cut -d'_' -f 1)
    setA["$SUBSTRING"]=1
done
printf '%s\n' "${!setA[@]}"

for m in "${!setA[@]}"; do 
mkdir -p "$folderpath/$m"
mv "$folderpath/${m}"* "$folderpath/${m}"
done




