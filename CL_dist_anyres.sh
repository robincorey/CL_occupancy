#!/bin/bash
# Script to measure dist between key residues and CL head groups
# Modified from time spent for CL sims 
## REQUIRES MANUAL TRIMMING OF XTC TO MATCH BILAYER FORMATION!!

ARRAY=(CG_0 CG_1 CG_2 CG_3 CG_4 CG_5 CG_6 CG_7 CG_8 CG_9)
ARRAYLEN=`echo ${#ARRAY[@]} - 1 | bc`
CURRENTDIR=`pwd`
RES=$1

if [ -z "${RES}" ]; then
        echo "add a fcuking variable you tard"
        exit 0
fi

####### Rest should be set in stone

for j in `seq 0 ${ARRAYLEN}`; do

cd $CURRENTDIR
cd ${ARRAY[j]} 

rm -f equilibration_3_${ARRAY[j]}.input.gro

TPR=equilibration_3_${ARRAY[j]}.tpr
#$GMX5/editconf -f equilibration_2.gro -o equilibration_3_${ARRAY[j]}.input.gro >& edc
GRO=equilibration_2.gro

rm -f index.ndx

echo -e "aPO1 | aPO2 | aGL0" '\n' q | make_ndx -f $GRO -o index.ndx >& md1

## First, define chains
chain1=`grep " 1MET     BB" $GRO | awk '{print $3}'`
chain2=`grep " 1ALA     BB" $GRO | awk '{print $3}'`; chain1b=`echo $chain2 - 1 | bc`
chain3=`grep " 1GLU     BB" $GRO | awk '{print $3}'`; chain2b=`echo $chain3 - 1 | bc`
chain4=`grep " 1HIS     BB" $GRO | awk '{print $3}'`; chain3b=`echo $chain4 - 1 | bc`
chain4b=`grep " 817DSPE   NH3" $GRO | awk '{print $3}'`
chain4c=`echo "$chain4b - 1" | bc`

rm -f chains.ndx

echo -e a$chain1-$chain1b'\n'a$chain2-$chain2b'\n'a$chain3-$chain3b'\n'a$chain4-$chain4c'\n'q | make_ndx -f $GRO -o chains.ndx -n index.ndx >& md2

if [[ ! -f chains.ndx ]]; then echo index making failed; exit 0 ; fi

indices=`grep -c "\[" chains.ndx`
## Groups for each chain
PO4=`echo "$indices - 5" | bc`

rm -f timespent_${ARRAY[j]}_${RES}.xvg
#echo timespent > timespent_${ARRAY[j]}_${RES}.xvg

echo "starting distances for ${ARRAY[j]}..."

CHAIN=(A Y E G)
IND=(4 3 2 1)
CHL=`echo "${#CHAIN[@]} -1" | bc`

if [[ ! -d xvgs/ ]]; then
        mkdir xvgs/
fi

for k in `seq 0 $CHL`; do
        LETTER=`echo "$indices - ${IND[k]}" | bc`
        echo $LETTER | editconf -f $GRO -o ${CHAIN[k]}.pdb -n chains.ndx >& eda
	CHARRAY=(`grep "BB  ${RES}" ${CHAIN[k]}.pdb | awk '{print $5}'`)
	CHAL=`echo ${#CHARRAY[@]} - 1 | bc` #
	for i in `seq 0 $CHAL`; do
		echo  $RES ${CHARRAY[i]}
		rm -f temp_${CHAIN[k]}.ndx
		echo -e "$LETTER & r${CHARRAY[i]}"'\n'q | make_ndx -f $GRO -o temp_${CHAIN[k]} -n chains.ndx >& mdA
		if [[ ! -f temp_${CHAIN[k]}.ndx ]]; then
			echo "#################
			make_ndx failed!
			##############"
			exit 0
		fi
		echo -e $indices'\n'PO1_PO2_GL0 | $GMX5/g_mindist -f BILAYER_${ARRAY[j]}.xtc -s $TPR -n temp_${CHAIN[k]} -od ${CHAIN[k]}_${CHARRAY[i]}_${ARRAY[j]}_${RES}.xvg -nice -19 -tu us >& dist
		echo Sec${CHAIN[k]}_${CHARRAY[i]} >> timespent_${ARRAY[j]}_${RES}.xvg 
        	num=`egrep -v "\#|\@" ${CHAIN[k]}_${CHARRAY[i]}_${ARRAY[j]}_${RES}.xvg  | awk '{print $2}' | egrep "\-01" | egrep "^1|^2|^3|^4|^5|^6|^7" | wc -l`
        	tot=`egrep -v "\#|\@" ${CHAIN[k]}_${CHARRAY[i]}_${ARRAY[j]}_${RES}.xvg  | wc -l`
        	occ=`echo "scale=3 ; $num / $tot * 100" | bc`
        	echo $occ >> timespent_${ARRAY[j]}_${RES}.xvg
		mv ${CHAIN[k]}_${CHARRAY[i]}*_${RES}.xvg xvgs/
	done
done


done

if [[ ! -d Backups/ ]]; then
	mkdir Backups/
fi

mv *#* Backups/

grep Sec timespent_${ARRAY[j]}_${RES}.xvg > Sec_${RES}.xvg
paste -d ',' timespent_*_${RES}.xvg > timespentall_${RES}.xvg
egrep -v Sec timespentall_${RES}.xvg > values_${RES}.xvg
paste -d ',' Sec_${RES}.xvg values_${RES}.xvg > timespentall_${RES}.xvg
