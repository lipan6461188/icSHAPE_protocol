#!/bin/bash

###########################################
#
# A Simulation of icSHAPE Data Analysis Pipeline
#
###########################################


###########################################
##	Help Information
info="""
Usgae: $0 [clear|run] \n
- clear: clear my workstation \n
- run: simulation icshape Pipeline \n\n
There are 9 Steps in icSHAPE pipeline: \n
\t1. Reads Collapses \n
\t2. Reads Trimming \n
\t3. Bowtie Map \n
\t\t3.1 Bowtie build index \n
\t\t3.2 Reads Map \n
\t4. Calculate RPKM \n
\t5. Count RT/BD \n
\t6. Combine RT replicates \n
\t7. Normalize RT \n
\t8. calcErichment \n
\t9. Filter Erichment\n
"""
if [ -z $1 ];then
	echo -e $info
	exit
fi

###########################################

###########################################
##	Clear Work Station
if [ $1 = 'clear' ]; then
	if [ -e collapes ]; then rm -rf collapes; fi
	if [ -e trimming ]; then rm -rf trimming; fi
	if [ -e bowtie ]; then rm -rf bowtie; fi
	if [ -e 45SrRNA.icSHAPE ]; then rm -rf 45SrRNA.icSHAPE; fi
	if [ -e enrichment ]; then rm -rf enrichment; fi
	if [ -e index ]; then rm -rf index; fi
	if [ -e normalize ]; then rm -rf normalize; fi
	if [ -e rpkm ]; then rm -rf rpkm; fi
	if [ -e rt ]; then rm -rf rt; fi
	if [ -e trimming ]; then rm -rf trimming; fi
	if [ -e combineRT ]; then rm -rf combineRT; fi
	exit
fi
###########################################

###########################################
##	Simulation icSHAPE Pipeline
if [ $1 = 'run' ]; then
	BIN=./icSHAPE/scripts
	ADAPTER=~/icSHAPE/data/adapter/TruSeq2-SE.fa  
	D1=./DATA/fq/D1.fq
	D2=./DATA/fq/D2.fq
	N1=./DATA/fq/N1.fq
	N2=./DATA/fq/N2.fq
	FA=./DATA/fa/45SrRNA.fa


	##	1. Reads Collapses
	if [ -e collapes ]; then rm -rf collapes; fi
	mkdir collapes

	$BIN/readCollapse.pl -U $D1 -o ./collapes/D1.fq
	$BIN/readCollapse.pl -U $D2 -o ./collapes/D2.fq
	$BIN/readCollapse.pl -U $N1 -o ./collapes/N1.fq
	$BIN/readCollapse.pl -U $N2 -o ./collapes/N2.fq


	##	2. Reads Trimming
	if [ -e trimming ]; then rm -rf trimming; fi
	mkdir trimming
	$BIN/trimming.pl -U ./collapes/D1.fq -o ./trimming/D1.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
	$BIN/trimming.pl -U ./collapes/D2.fq -o ./trimming/D2.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
	$BIN/trimming.pl -U ./collapes/N1.fq -o ./trimming/N1.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
	$BIN/trimming.pl -U ./collapes/N2.fq -o ./trimming/N2.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25


	##	3. Bowtie Map

	###		3.1 Bowtie build index 
	if [ -e index ]; then rm -rf index; fi
	mkdir index
	bowtie2-build $FA index/rRNA

	###		3.2 Reads Map
	if [ -e bowtie ]; then rm -rf bowtie; fi
	mkdir bowtie
	bowtie2 -U ./trimming/D1.fq -S ./bowtie/D1.sam -x ./index/rRNA --non-deterministic -p 2
	bowtie2 -U ./trimming/D2.fq -S ./bowtie/D2.sam -x ./index/rRNA --non-deterministic -p 2
	bowtie2 -U ./trimming/N1.fq -S ./bowtie/N1.sam -x ./index/rRNA --non-deterministic -p 2
	bowtie2 -U ./trimming/N2.fq -S ./bowtie/N2.sam -x ./index/rRNA --non-deterministic -p 2


	##	4. Calculate RPKM
	if [ -e rpkm ]; then rm -rf rpkm; fi
	mkdir rpkm
	$BIN/estimateRPKM.pl -i ./bowtie/D1.sam -o ./rpkm/D1.rpkm -c 2
	$BIN/estimateRPKM.pl -i ./bowtie/D2.sam -o ./rpkm/D2.rpkm -c 2
	$BIN/estimateRPKM.pl -i ./bowtie/N1.sam -o ./rpkm/N1.rpkm -c 2
	$BIN/estimateRPKM.pl -i ./bowtie/N2.sam -o ./rpkm/N2.rpkm -c 2


	##	5. Count RT/BD
	if [ -e rt ]; then rm -rf rt; fi
	mkdir rt
	$BIN/calcRT.pl -i ./bowtie/D1.sam -o ./rt/D1.rt -r ./rpkm/D1.rpkm -c 2
	$BIN/calcRT.pl -i ./bowtie/D2.sam -o ./rt/D2.rt -r ./rpkm/D2.rpkm -c 2
	$BIN/calcRT.pl -i ./bowtie/N1.sam -o ./rt/N1.rt -r ./rpkm/N1.rpkm -c 2
	$BIN/calcRT.pl -i ./bowtie/N2.sam -o ./rt/N2.rt -r ./rpkm/N2.rpkm -c 2

	##	6. Combine RT replicates
	if [ -e combineRT ]; then rm -rf combineRT; fi
	mkdir combineRT
	$BIN/combineRTreplicates.pl -i ./rt/D1.rt:./rt/D2.rt -o ./combineRT/D.rt
	$BIN/combineRTreplicates.pl -i ./rt/N1.rt:./rt/N2.rt -o ./combineRT/N.rt

	##	7. Normalize RT
	if [ -e normalize ]; then rm -rf normalize; fi
	mkdir normalize
	$BIN/normalizeRTfile.pl -i ./combineRT/D.rt -o ./normalize/D.norm.rt -m mean:vigintile2 -d 32 -l 32 -f 100
	$BIN/normalizeRTfile.pl -i ./combineRT/N.rt -o ./normalize/N.norm.rt -m mean:vigintile2 -d 32 -l 32 -f 100

	##	8. calcErichment
	if [ -e enrichment ]; then rm -rf enrichment; fi
	mkdir enrichment
	$BIN/calcEnrich.pl -f ./normalize/N.norm.rt -b ./normalize/D.norm.rt -o ./enrichment/45SrRNA.erichment -w factor5:scaling1 -y 10 -x 0.25 -e complex

	##	9. Filter Erichment
	$BIN/filterEnrich.pl -i ./enrichment/45SrRNA.erichment -o 45SrRNA.icSHAPE -t 100 -T 2 -s 5 -e 30
fi
###########################################


