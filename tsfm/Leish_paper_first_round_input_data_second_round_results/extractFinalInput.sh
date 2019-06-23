#!/bin/bash

folderpath=$(pwd) 
# list all the folders except for trash
folders=$(ls -l $folderpath | grep "^d" | awk -F" " '{print $9}')

for name in $folders; do 
if [ $name != "trash" ] && [ $name != "first_round_results" ] && [ $name != "HOMO" ] && [ $name != "LogotaxAlignmentInputs" ]
then
mkdir -p "$folderpath/LogotaxAlignmentInputs/sites74/$name/"
mkdir -p "$folderpath/LogotaxAlignmentInputs/sites72/$name/"
mkdir -p "$folderpath/LogotaxAlignmentInputs/KLD/$name/"
mkdir -p "$folderpath/LogotaxAlignmentInputs/ID/$name/"
mkdir -p "$folderpath/LogotaxAlignmentInputs/Func/$name/"

mv "$folderpath/$name/"*sites74*.aln "$folderpath/LogotaxAlignmentInputs/sites74/$name/"
mv "$folderpath/$name/"*sites72*.aln "$folderpath/LogotaxAlignmentInputs/sites72/$name/"
mv "$folderpath/$name/"*KLD*.eps "$folderpath/LogotaxAlignmentInputs/KLD/$name/"
mv "$folderpath/$name/"*ID*.eps "$folderpath/LogotaxAlignmentInputs/ID/$name/"
mv "$folderpath/$name/"*.eps "$folderpath/LogotaxAlignmentInputs/Func/$name/"

fi
done

