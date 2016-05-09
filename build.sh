#!/bin/bash

if [ -z $1 ]  
then
	run=20 
else
	run=$1
fi

echo "$run runs per model"

#Number of bottlenecks
bottleneckcount="50"
#Number of generations
growthtime="25"
#Growth limit when the simulation exceeds this barier bottlenecks to 1
limgrowth="1000000"

#Total number of genes
numberofgenes="4400"
#Genome size
genomesize="4640000"
#number of genes which may change the mutation rates
numberofmutatorloci="20"

#mutation rate for deleterious loci
mutationrate_deleterious="0.9"
#mutation rate for beneficial loci
mutationrate_beneficial="0.01"
#mutation rate for neutral loci
mutationrate_neutral="9.9"
#at each bottleneck an anti mutator appears with this probability
#increasedmutatormutationrate="0.5"
#the multiplication factor which changes the mutation rate
meaneffectmutator="0.9"
varianceeffectmutator="0.2"
#the constants of fitness to offspring conversion equation
wtooff0="-0.26"
wtooff1="8.54"
wtooff2="-10.26"
wtooff3="4.34"
#the amount of variation in the offspring number due to the fitness function
fitnessselection="0.25"

#deleterious mutations' fitness effects are distributed with a gamma distribution shape and scale factors
shapeofgammadel="0.3"
scaleofgammadel="0.1"
#beneficial mutations' fitness effects are distributed with an exponential distribution the beta value of it
betaben="0.01"

#the highest percentage of the population which could be selected as the parent of next bottleneck
selecthighest="0.5"


echo "#ifndef __VARIABLES_H_" 					>> variables.tmp
echo "#define __VARIABLES_H_" 					>> variables.tmp
echo " " 							>> variables.tmp
echo " " 							>> variables.tmp
echo "#define T_bot            $bottleneckcount" 		>> variables.tmp
echo "#define T_growth         $growthtime" 			>> variables.tmp
echo "#define N_growth         $limgrowth"			>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Ngenes           $numberofgenes" 			>> variables.tmp
echo "#define L                $genomesize" 			>> variables.tmp
echo "#define Nmutator         $numberofmutatorloci"		>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Mu_d             $mutationrate_deleterious"	>> variables.tmp
echo "#define Mu_b             $mutationrate_beneficial" 	>> variables.tmp
echo "#define Mu_n             $mutationrate_neutral" 		>> variables.tmp
#echo "#define Inc_mut_rate     $increasedmutatormutationrate" 	>> variables.tmp
echo "#define Mean_e_mutator   $meaneffectmutator" 		>> variables.tmp
echo "#define Var_e_mutator    $varianceeffectmutator"		>> variables.tmp
echo "#define WtoOffa          $wtooff0"			>> variables.tmp
echo "#define WtoOffb          $wtooff1"			>> variables.tmp
echo "#define WtoOffc          $wtooff2"			>> variables.tmp
echo "#define WtoOffd          $wtooff3"			>> variables.tmp
echo "#define W_sel            $fitnessselection"               >> variables.tmp
echo " " 							>> variables.tmp
echo "#define Shaped           $shapeofgammadel" 		>> variables.tmp
echo "#define Scaled           $scaleofgammadel" 		>> variables.tmp
echo "#define Betab            $betaben"	 		>> variables.tmp
echo " " 							>> variables.tmp
echo "#define Selecthighest    $selecthighest" 			>> variables.tmp
echo " " 							>> variables.tmp

mkdir	sims

mkdir 	sims/nomuttr
cp 	variables.tmp	variables.h
echo "#define Mutator_control      0" 	>> variables.h
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

rm variables.*
