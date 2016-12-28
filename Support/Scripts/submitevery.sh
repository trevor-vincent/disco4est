#!/bin/bash 

NOIDCOUNTER=0
IDCOUNTER=0
for i in $( ls -d */ ); do  
    cd $i 
    if [ ! -e "job.id" ] 
    then 
	((NOIDCOUNTER++))
    else
	((IDCOUNTER++))
    fi 
    cd ..
done
echo "Folders with jobid = $IDCOUNTER, Folders without jobid = $NOIDCOUNTER"

STOPNUM=100
USER=tvincent
COUNTER=0
for i in $( ls -d */ ); do  
    cd $i 
    if [ ! -e "job.id" ] 
    then 
	((COUNTER++))
	echo $COUNTER
	qstat -u ${USER} > temp0.txt 
	qsub submit.sh 
	qstat -u ${USER} > temp1.txt
	diff temp0.txt temp1.txt > job.id
	rm -f temp0.txt
	rm -f temp1.txt	
    fi 
    if [ $COUNTER -gt $STOPNUM ]
    then
	exit
    fi
    cd ..
done

$NOIDCOUNTER2=0
$IDCOUNTER2=0
for i in $( ls -d */ ); do  
    cd $i 
    if [ ! -e "job.id" ] 
    then 
	((NOIDCOUNTER2++))
    else
	((IDCOUNTER2++))
    fi 
    cd ..
done
echo "Folders with jobid = $IDCOUNTER2, Folders without jobid = $NOIDCOUNTER2"
