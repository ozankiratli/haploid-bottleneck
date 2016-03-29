#!/bin/bash

if [ -z $1 ]  
then
	run=20 
else
	run=$1
fi

echo "$run runs per model"

bottleneckcount="50"
growthtime="25"

numberofgenes="4400"
genomesize="4640000"
numberofmutatorloci="20"

mutationrate_deleterious="0.9"
mutationrate_beneficial="0.001"
mutationrate_neutral="9.9"
increasedmutatormutationrate="0.1"
meaneffectofmutator="0.9"

shapeofgammadel="0.3"
scaleofgammadel="0.1"
betaben="0.1"

selecthighest="0.5"


echo "#ifndef __VARIABLES_H_" 					>> variables.tmp
echo "#define __VARIABLES_H_" 					>> variables.tmp
echo " " 							>> variables.tmp
echo " " 							>> variables.tmp
echo "#define T_bot            $bottleneckcount" 		>> variables.tmp
echo "#define T_growth         $growthtime" 			>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Ngenes           $numberofgenes" 			>> variables.tmp
echo "#define L                $genomesize" 			>> variables.tmp
echo "#define Nmutator         $numberofmutatorloci"		>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Mu_d             $mutationrate_deleterious"	>> variables.tmp
echo "#define Mu_b             $mutationrate_beneficial" 	>> variables.tmp
echo "#define Mu_n             $mutationrate_neutral" 		>> variables.tmp
echo "#define Inc_mut_rate     $increasedmutatormutationrate" 	>> variables.tmp
echo "#define Mean_e_mutator   $meaneffectofmutator" 		>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Shaped           $shapeofgammadel" 		>> variables.tmp
echo "#define Scaled           $scaleofgammadel" 		>> variables.tmp
echo "#define Betab            $betaben"	 		>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Selecthighest            $selecthighest" 			>> variables.tmp
echo " " 							>> variables.tmp

mkdir	sims

mkdir 	sims/nomuttr
cp 	variables.tmp	variables.h
echo "#define Mutator_control      0" 	>> variables.h
echo "#define Increased_mutator    0" 	>> variables.h	
echo "#define Del_dep              0"	>> variables.h
echo " " 				>> variables.h

for ((runs=0;runs<$run;runs++))
do
	mkdir	sims/nomuttr/run$runs
	cp	simulation.cpp	sims/nomuttr/run$runs
	cp	variables.h	sims/nomuttr/run$runs
	echo "#define Run         $runs" 		>> sims/nomuttr/run$runs/variables.h
	echo " "    					>> sims/nomuttr/run$runs/variables.h
	echo "#endif"					>> sims/nomuttr/run$runs/variables.h
	g++ -std=c++11 sims/nomuttr/run$runs/simulation.cpp -o sims/nomuttr/run$runs/simulation
done



mkdir sims/mutator
cp variables.tmp	variables.h
echo "#define Mutator_control      1" 	>> variables.h
echo "#define Increased_mutator    0" 	>> variables.h
echo "#define Del_dep              0"	>> variables.h
echo " " 				>> variables.h
for ((runs=0;runs<$run;runs++))
do
	mkdir   sims/mutator/run$runs
	cp      simulation.cpp   sims/mutator/run$runs
	cp      variables.h    sims/mutator/run$runs
        echo "#define Run               $runs"		>> sims/mutator/run$runs/variables.h
       	echo " "    					>> sims/mutator/run$runs/variables.h
	echo "#endif"					>> sims/mutator/run$runs/variables.h
	g++ -std=c++11 sims/mutator/run$runs/simulation.cpp -o sims/mutator/run$runs/simulation
done



mkdir   sims/incmttr
cp      variables.tmp   variables.h
echo "#define Mutator_control      1" 	>> variables.h
echo "#define Increased_mutator    1" 	>> variables.h	
echo "#define Del_dep              0"	>> variables.h
echo " " 				>> variables.h
for ((runs=0;runs<$run;runs++))
do
	mkdir   sims/incmttr/run$runs
	cp      simulation.cpp   sims/incmttr/run$runs
	cp      variables.h      sims/incmttr/run$runs
        echo "#define Run               $runs"		>> sims/incmttr/run$runs/variables.h
       	echo " "    					>> sims/incmttr/run$runs/variables.h
	echo "#endif"					>> sims/incmttr/run$runs/variables.h
	g++ -std=c++11 sims/incmttr/run$runs/simulation.cpp -o sims/incmttr/run$runs/simulation
done



mkdir   sims/deldepn
cp      variables.tmp   variables.h
echo "#define Mutator_control      1" 	>> variables.h
echo "#define Increased_mutator    1" 	>> variables.h	
echo "#define Del_dep              1"	>> variables.h
echo " " 				>> variables.h
for ((runs=0;runs<$run;runs++))
do
	mkdir   sims/deldepn/run$runs
	cp      simulation.cpp   sims/deldepn/run$runs
	cp      variables.h      sims/deldepn/run$runs
	echo "#define Run               $runs"		>> sims/deldepn/run$runs/variables.h
	echo " "    					>> sims/deldepn/run$runs/variables.h
	echo "#endif"					>> sims/deldepn/run$runs/variables.h
	g++ -std=c++11 sims/deldepn/run$runs/simulation.cpp -o sims/deldepn/run$runs/simulation
done

rm variables.*
