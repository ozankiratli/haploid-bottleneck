#!/bin/bash

if [ "$1" == "" ]  
then
	x=1
else
	x=$1
fi

echo "$x parallel simulations at a time"


for dirs1 in  $PWD/sims/* ; do
	cd $dirs1
	for dirs2 in $PWD/* ;
	do
		echo $dirs2
		cd $dirs2
		nohup ./simulation  &> simulation.log&
		sleep 10s
		while [ `ps aux -U ozan | grep simulation | wc -l` -gt "$x" ] 	
		do
			echo "Sleeping for 1min until simulation count drops under $x"
			sleep 1m
		done
		cd ..
	done
	cd ..
done
