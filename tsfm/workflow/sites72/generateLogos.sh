#!/bin/bash
# this script will run tsfm on all the classes of TriTryp genomes
# for each genomes it will make a folder in each folder will make three folders: KLD,ID,Func, bubble_Table


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
	mkdir -p "$folderpath/Logos/$name1/Func"
	mkdir -p "$folderpath/Logos/$name1/Bubble"
#run tsfm to make the logos 
name2=${name1%".v5"}
name3="$name2.sites72.v5"
#printf "$name3"
python3 "$tsfmpath/tsfm.py" -c "$folderpath/tRNA_L_skel_Leish_sites72_struct.txt" --logo "$folderpath/$name1/$name3"
mv -- *.eps "$folderpath/Logos/$name1/Func"
fi
done

# creating ID, KLD, and tables for the bubble plots
for name1 in $folders; do 
if [ $name1 != "Logos" ]
then
printf "***************$name1***************"
name2=${name1%".v5"}
name3="$name2.sites72.v5"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/tRNA_L_skel_Leish_sites72_struct.txt" --idlogo "$folderpath/$name1/$name3" "$folderpath/HOMO.v5/HOMO.sites72.v5"
mv -- *.eps "$folderpath/Logos/$name1/ID"

python3 "$tsfmpath/tsfm.py" -c "$folderpath/tRNA_L_skel_Leish_sites72_struct.txt" --kldlogo --bt "$folderpath/$name1/$name3" "$folderpath/HOMO.v5/HOMO.sites72.v5"
mv -- *.eps "$folderpath/Logos/$name1/KLD"
mv *_Table.txt "$folderpath/Logos/$name1/Bubble"

fi
done



