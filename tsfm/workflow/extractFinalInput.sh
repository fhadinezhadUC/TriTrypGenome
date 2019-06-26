#!/bin/bash

folderpath=$(pwd) 
# list all the folders except for trash
folders=$(ls -l $folderpath | grep "^d" | awk -F" " '{print $9}')

for name in $folders; do 
if [ $name != "trash" ] && [ $name != "first_round_results" ] && [ $name != "HOMO" ]
then
mkdir -p "$folderpath/LogotaxAlignmentInputs/sites74/$name/"
mkdir -p "$folderpath/LogotaxAlignmentInputs/sites72/$name/"
mv "$folderpath/$name/"*sites74*.aln "$folderpath/LogotaxAlignmentInputs/sites74/$name/"
mv "$folderpath/$name/"*sites72*.aln "$folderpath/LogotaxAlignmentInputs/sites72/$name/"
fi
done

