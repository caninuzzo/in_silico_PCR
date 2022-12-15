#!/bin/bash

while getopts ":s:p:" opt
   do
     case $opt in
        s ) seq=$OPTARG;;
        p ) primers=$OPTARG;;

     esac
done

primersOnly=$(echo "$primers" | rev | cut -d"/" -f1 | rev )

F=$(echo "$primersOnly" | cut -d'_' -f1)
R=$(echo "$primersOnly" | cut -d'_' -f2)

echo "_______________________________"
echo -e "\n Amorces choisies : \n Forward : $F // Reverse : $R \n"
echo "_______________________________"

### Choisir le nombre maximum de mismatches autorisé

echo -e "\n Choisir le nombre maximum de mismatches autorisé : \n"
echo -e "\n ------------------------------- \n"
read diffs


### Ouverture de mothur 

mothur "#pcr.seqs(fasta=$seq, oligos=$primers, pdiffs=$diffs, rdiffs=$diffs, keepdots=F, keepprimer=T)"

### for displaying results : 
nb_seq=$(grep -c "^>" $seq)

amplicon_file=$(echo $seq | sed "s/.fasta$/.pcr.fasta/g")
nb_seq_ampl=$(grep -c "^>" $amplicon_file)

### get summary of the amplicons :

mothur "#summary.seqs(fasta=$amplicon_file)" 1> fmr_sum_seq

grep -A10 "Start" fmr_sum_seq

let Median=$(grep "Median" fmr_sum_seq | cut -f4)

let CutMin=$(echo $Median*0.8 | bc | awk '{print int($0)}')

let CutMax=$(echo $Median*1.2 | bc | awk '{print int($0)}')

echo "Median = $Median / CutMin = $CutMin / CutMax = $CutMax"

echo -e "\n$nb_seq_ampl/$nb_seq séquences ont été sélectionnées grâce au couple d'amorce ($F/$R) avec un taux d'erreur max autorisé de $diffs\n Voici un résumé de celles-ci : \n"

#Création de la matrice pour les amorces

# find motifs in pcr.fasta to adapt the analysis :
Ftest=$(grep -c "fpdiffs" $amplicon_file)
Rtest=$(grep -c "rpdiffs" $amplicon_file)

### If we test only reverse primer :
if [[ $Ftest -eq 0 ]] && [[ $Rtest -ne 0 ]];#analyse reverse only 
then
declare -a matrix
matrix[0]="rdiffs"
#count:
vec[0]="Nb seq"
for i in `seq 0 $diffs` ;
do
matrixWrite=";x=$i"
matrix[0]="${matrix[0]} ${matrixWrite}"
##
vec[$i+1]=$(grep -c "rpdiffs=$i" $amplicon_file)
done
#assemble:
matrix[1]=$(IFS=";"; echo "${vec[*]}")
fi

### If we test only forward primer :
if [[ $Ftest -ne 0 ]] && [[ $Rtest -eq 0 ]];#analyse forward only 
then
echo "FORWARD ONLY !"
declare -a matrix
matrix[0]="fdiffs"
#count:
vec[0]="Nb seq"
for i in `seq 0 $diffs` ;
do
matrixWrite=";x=$i"
matrix[0]="${matrix[0]} ${matrixWrite}"
##
vec[$i+1]=$(grep -c "fpdiffs=$i" $amplicon_file)
done
#assemble:
matrix[1]=$(IFS=";"; echo "${vec[*]}")

outMatrix=${matrix[0]}
outMatrix="${outMatrix}\n${matrix[1]}"

fi

### If we test a primers couple :
if [[ $Ftest -ne 0 ]] && [[ $Rtest -ne 0 ]];#analyse both primers
then
declare -a matrix
matrix[0]="fdiffs vs rdiffs"
# double counts:
for i in `seq 0 $diffs` ;
do
matrixWrite=";x=$i"
matrix[0]="${matrix[0]} ${matrixWrite}"
vec[0]="x=$i"
for j in `seq 0 $diffs` ;
do
vec[$j+1]=$(grep -c "fpdiffs=$i(match) rpdiffs=$j" $amplicon_file)
#assemble
if [ $j -eq $diffs ]
then
matrix[$i+1]=$(IFS=";"; echo "${vec[*]}")
fi
done #j
done #i

### header :
outMatrix=${matrix[0]}
for i in `seq 0 $diffs` ;
do
outMatrix="${outMatrix}\n${matrix[$i+1]}"
done
fi

### If we've tested something but nothing matched :
if [[ $Ftest -eq 0 ]] && [[ $Rtest -eq 0 ]];#nothing has been found -> exit
then
echo -e "\n No sequences seems to match with this primer or primer couple. Try with Other\n Bye........... !\n "
exit
fi


echo -e "Voici un aperçu de la matrice créée : \n$outMatrix\n"

echo -e "\n .............................. Bye.\n"

####### suite : spécific resolution :

mothur "#screen.seqs(fasta=$amplicon_file,minlength=$CutMin,maxlength=$CutMax)" 1> fmr.screen

amplicon_screened_file=$(echo $seq | sed "s/.fasta$/.pcr.good.fasta/g")

totAmp=$(grep -c "^>" $amplicon_screened_file)

echo "nb total amplifié : $totAmp"

mothur "#unique.seqs(fasta=$amplicon_screened_file)" 1> fmr.uniq

amplicon_derep_file=$(echo $seq | sed "s/.fasta$/.pcr.good.unique.fasta/g")

Derep=$(grep -c "^>" $amplicon_derep_file)

echo "amplicon uniques récupérés : $Derep"


Resol=$(echo "scale=4;$Derep/$totAmp*100" | bc )

echo -e "\n Résolution de cette région : $Resol % \n"


#on suppprime les fichiers mtn inutiles
pathFile=$(echo $seq | rev | cut -d"/" -f2- | rev )

rm $pathFile/*.bad.accnos
rm $pathFile/*.scrap.pcr.fasta
rm *.logfile
rm fmr*