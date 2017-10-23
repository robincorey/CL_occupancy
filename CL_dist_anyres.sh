#!/bin/bash
###############################################################
# Script to measure dist between all residues in a protein MD #
# trajectory and a lipid head group (can define). MARTINI     #
###############################################################

TRAJ=$1    ; TRAJ2=${TRAJ::-4}
TPR=$2     ; TPR2=${TPR::-4}
RES=$3   
LIPID=$4

if [ "${TRAJ}" == "-h" ]; then
	echo "
Run this programme in your working folder, and supply input files as follows:
./CL_dist_anyres.sh trajectoryname.xtc runinput.tpr LYS
"
	exit 0
fi

if [ -z "${TRAJ}" ]; then
        echo "Need to supply a trajectory name
run ./CL_dist_anyres.sh -h for help"
        exit 0
fi

if [ ${TRAJ: -3} != "xtc" ] && [ ${TRAJ: -3} != "trr" ]; then
        echo "Trajectory needs to be in trr or xtc format
run ./CL_dist_anyres.sh -h for help"
        exit 0
fi

if [ ${TPR: -3} != "tpr" ]; then
        echo "Tpr file needed
run ./CL_dist_anyres.sh -h for help"
        exit 0
fi

if [ -z "${RES}" ]; then
        echo "Please supply a three character residue name, e.g. LYS
run ./CL_dist_anyres.sh -h for help"
        exit 0
fi

rm -f ${TPR2}.input.gro
editconf -f ${TPR} -o ${TPR2}.input.gro >& edc_output
GRO=${TPR2}.input.gro

rm -f index.ndx
phos=($(grep CAR ${TPR2}.input.gro  | grep " P"  | awk '{print $3}' | sort -u))
echo -e "aPO1 | aPO2 | aGL0" '\n' q | make_ndx -f $GRO -o index.ndx >& md1

## First, define chains. Currently 3DIN specific
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
echo timespent > timespent_${ARRAY[j]}_${RES}.xvg

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
