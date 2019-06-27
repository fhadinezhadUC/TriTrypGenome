#!/bin/bash
# this script will run tsfm on all the classes of TriTryp genomes
# for each genomes it will make a folder in each folder will make three folders: KLD,ID,bubble_Table


folderpath=$(pwd) 
tsfmpath="/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfm-master/tsfm"
folders=$(ls -l $folderpath | grep "^d" | awk -F" " '{print $9}')

# creating function logos
for name1 in $folders; do 
if [ $name1 != "Logos" ]
then
printf "***************$name1***************"
	mkdir -p "$folderpath/Logos/$name1/KLD"
	mkdir -p "$folderpath/Logos/$name1/ID"
	#mkdir -p "$folderpath/Logos/$name1/Func"
	mkdir -p "$folderpath/Logos/$name1/Bubble"
python3 "$tsfmpath/tsfm.py" -c "$folderpath/tRNA_L_skel_Leish.sites74.struct.txt" --idlogo --kldlogo --bt "$folderpath/$name1/$name1" "$folderpath/HOMO/HOMO"
mv -- *ID*.eps "$folderpath/Logos/$name1/ID"
mv -- *KLD*.eps "$folderpath/Logos/$name1/KLD"
mv *_Table.txt "$folderpath/Logos/$name1/Bubble"
fi
done




