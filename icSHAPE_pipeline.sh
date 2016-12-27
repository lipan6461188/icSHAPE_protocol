#!/bin/bash

###########################################
#
# A Simulation of icSHAPE Data Analysis Pipeline
#
###########################################


###########################################
##	Help Information
info="""
Usgae: $0 [clear|run|cov|step(n)] \n
- clear: clear my workstation \n
- run: simulate icshape Pipeline \n
- cov: converse icshape to VARNA color map file format \n
- step: select a step below to run \n
\t$0 step1 [D1.fq D2.fq N1.fq N2.fq]\n
\t$0 step2 [D1.collaped.fq D2.collaped.fq N1.collaped.fq N2.collaped.fq]\n
\t$0 step3 [D1.trimed.fq D2.trimed.fq N1.trimed.fq N2.trimed.fq fasta]\n
\t$0 step4 [D1.sam D2.sam N1.sam N2.sam]\n
\t$0 step5 [D1.sam D2.sam N1.sam N2.sam D1.rpkm D2.rpkm N1.rpkm N2.rpkm]\n
\t$0 step6 [D1.rt D2.rt N1.rt N2.rt]\n
\t$0 step7 [D.rt N.rt]\n
\t$0 step8 [D.norm.rt N.norm.rt]\n
\t$0 step9 [rRNA.erichment rRNA.icSHAPE]\n
     Example: $0 step5 bowtie/D1.sam bowtie/D2.sam bowtie/N1.sam bowtie/N2.sam rpkm/D1.rpkm rpkm/D2.rpkm rpkm/N1.rpkm rpkm/N2.rpkm\n\n
There are 9 Steps in icSHAPE pipeline: \n
\t[1]. Reads Collapses \n
\t[2]. Reads Trimming \n
\t[3]. Bowtie Map \n
\t\t[3.1] Bowtie build index \n
\t\t[3.2] Reads Map \n
\t[4]. Calculate RPKM \n
\t[5]. Count RT/BD \n
\t[6]. Combine RT replicates \n
\t[7]. Normalize RT \n
\t[8]. calcErichment \n
\t[9]. Filter Erichment\n
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
	if [ -e rRNA.icSHAPE ]; then rm -rf rRNA.icSHAPE; fi
	if [ -e enrichment ]; then rm -rf enrichment; fi
	if [ -e index ]; then rm -rf index; fi
	if [ -e normalize ]; then rm -rf normalize; fi
	if [ -e rpkm ]; then rm -rf rpkm; fi
	if [ -e rt ]; then rm -rf rt; fi
	if [ -e trimming ]; then rm -rf trimming; fi
	if [ -e combineRT ]; then rm -rf combineRT; fi
	if [ -e rRNA.icSHAPE.varna ]; then rm rRNA.icSHAPE.varna; fi
	exit
fi
###########################################

###########################################
##	Simulation icSHAPE Pipeline

#FA=./DATA/fa/5SrRNA.fa
#FA=./DATA/fa/18SrRNA.fa
#FA=./DATA/fa/45SrRNA.fa

if [ -z "$(uname -a | grep -i 'Linux')" ];
then
    echo 'OS: Mac OS';
    step_num=$(echo $1 | grep -E "^step\d+$" | sed 's/step//')
else
    echo 'OS: Linux'
    step_num=$(echo $1 | grep -P "^step\d+$" | sed 's/step//')
fi

if [ $1 = 'run' ] || [ ! -z $step_num ]; then
	BIN=./icSHAPE/scripts
	ADAPTER=./icSHAPE/data/adapter/TruSeq2-SE.fa  

	##	1. Reads Collapses
	if [ $1 = 'run' ] || [ $step_num = 1 ]; then
	    if [ -z $2 ]; then D1=./DATA/fq/D1.fq; else D1=$2; fi
		if [ -z $3 ]; then D2=./DATA/fq/D2.fq; else D2=$3; fi
		if [ -z $4 ]; then N1=./DATA/fq/N1.fq; else N1=$4; fi
		if [ -z $5 ]; then N2=./DATA/fq/N2.fq; else N2=$5; fi

		if [ -e collapes ]; then rm -rf collapes; fi
		mkdir collapes
		$BIN/readCollapse.pl -U $D1 -o ./collapes/D1.fq
		$BIN/readCollapse.pl -U $D2 -o ./collapes/D2.fq
		$BIN/readCollapse.pl -U $N1 -o ./collapes/N1.fq
		$BIN/readCollapse.pl -U $N2 -o ./collapes/N2.fq
	fi

	##	2. Reads Trimming
	if [ $1 = 'run' ] || [ $step_num = 2 ]; then
		if [ -z $2 ]; then D1=./collapes/D1.fq; else D1=$2; fi                        
        if [ -z $3 ]; then D2=./collapes/D2.fq; else D2=$3; fi                        
        if [ -z $4 ]; then N1=./collapes/N1.fq; else N1=$4; fi                        
        if [ -z $5 ]; then N2=./collapes/N2.fq; else N2=$5; fi                        
		if [ -z $6 ]; then ADAPTER=./icSHAPE/data/adapter/TruSeq2-SE.fa; else ADAPTER=$6; fi

		if [ -e trimming ]; then rm -rf trimming; fi
		mkdir trimming
		$BIN/trimming.pl -U $D1 -o ./trimming/D1.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
		$BIN/trimming.pl -U $D2 -o ./trimming/D2.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
		$BIN/trimming.pl -U $N1 -o ./trimming/N1.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
		$BIN/trimming.pl -U $N2 -o ./trimming/N2.fq -l 13 -t 0 -c phred33 -a $ADAPTER -m 25
	fi

	##	3. Bowtie Map
	if [ $1 = 'run' ] || [ $step_num = 3 ]; then
        if [ -z $2 ]; then D1=./trimming/D1.fq; else D1=$2; fi                         
        if [ -z $3 ]; then D2=./trimming/D2.fq; else D2=$3; fi                         
        if [ -z $4 ]; then N1=./trimming/N1.fq; else N1=$4; fi                         
        if [ -z $5 ]; then N2=./trimming/N2.fq; else N2=$5; fi 
		if [ -z $6 ]; then FA=./DATA/fa/18SrRNA.fa; else FA=$6; fi

		###		3.1 Bowtie build index 
		if [ -e index ]; then rm -rf index; fi
		mkdir index
		bowtie2-build $FA index/rRNA

		###		3.2 Reads Map
		if [ -e bowtie ]; then rm -rf bowtie; fi
		mkdir bowtie
		bowtie2 -U $D1 -S ./bowtie/D1.sam -x ./index/rRNA --non-deterministic -p 9
		bowtie2 -U $D2 -S ./bowtie/D2.sam -x ./index/rRNA --non-deterministic -p 9
		bowtie2 -U $N1 -S ./bowtie/N1.sam -x ./index/rRNA --non-deterministic -p 9
		bowtie2 -U $N2 -S ./bowtie/N2.sam -x ./index/rRNA --non-deterministic -p 9
	fi

	##	4. Calculate RPKM
	if [ $1 = 'run' ] || [ $step_num = 4 ]; then
        if [ -z $2 ]; then D1=./bowtie/D1.sam; else D1=$2; fi                     
        if [ -z $3 ]; then D2=./bowtie/D2.sam; else D2=$3; fi                     
        if [ -z $4 ]; then N1=./bowtie/N1.sam; else N1=$4; fi                     
        if [ -z $5 ]; then N2=./bowtie/N2.sam; else N2=$5; fi       

		if [ -e rpkm ]; then rm -rf rpkm; fi
		mkdir rpkm
		$BIN/estimateRPKM.pl -i $D1 -o ./rpkm/D1.rpkm -c 2
		$BIN/estimateRPKM.pl -i $D2 -o ./rpkm/D2.rpkm -c 2
		$BIN/estimateRPKM.pl -i $N1 -o ./rpkm/N1.rpkm -c 2
		$BIN/estimateRPKM.pl -i $N2 -o ./rpkm/N2.rpkm -c 2
	fi

	##	5. Count RT/BD
	if [ $1 = 'run' ] || [ $step_num = 5 ]; then
        if [ -z $2 ]; then D1=./bowtie/D1.sam; else D1=$2; fi
        if [ -z $3 ]; then D2=./bowtie/D2.sam; else D2=$3; fi
        if [ -z $4 ]; then N1=./bowtie/N1.sam; else N1=$4; fi
        if [ -z $5 ]; then N2=./bowtie/N2.sam; else N2=$5; fi

        if [ -z $6 ]; then rD1=./rpkm/D1.rpkm; else rD1=$6; fi
        if [ -z $7 ]; then rD2=./rpkm/D2.rpkm; else rD2=$7; fi
        if [ -z $8 ]; then rN1=./rpkm/N1.rpkm; else rN1=$8; fi
        if [ -z $9 ]; then rN2=./rpkm/N2.rpkm; else rN2=$9; fi


		if [ -e rt ]; then rm -rf rt; fi
		mkdir rt
		$BIN/calcRT.pl -i $D1 -o ./rt/D1.rt -r $rD1 -c 2
		$BIN/calcRT.pl -i $D2 -o ./rt/D2.rt -r $rD2 -c 2
		$BIN/calcRT.pl -i $N1 -o ./rt/N1.rt -r $rN1 -c 2
		$BIN/calcRT.pl -i $N2 -o ./rt/N2.rt -r $rN2 -c 2
	fi

	##	6. Combine RT replicates
        if [ -z $2 ]; then D1=./rt/D1.rt; else D1=$2; fi
        if [ -z $3 ]; then D2=./rt/D2.rt; else D2=$3; fi
        if [ -z $4 ]; then N1=./rt/N1.rt; else N1=$4; fi
        if [ -z $5 ]; then N2=./rt/N2.rt; else N2=$5; fi


	if [ $1 = 'run' ] || [ $step_num = 6 ]; then
		if [ -e combineRT ]; then rm -rf combineRT; fi
		mkdir combineRT
		$BIN/combineRTreplicates.pl -i $D1:$D2 -o ./combineRT/D.rt
		$BIN/combineRTreplicates.pl -i $N1:$N2 -o ./combineRT/N.rt
	fi

	##	7. Normalize RT
        if [ -z $2 ]; then D=./combineRT/D.rt; else D=$2; fi
        if [ -z $3 ]; then N=./combineRT/N.rt; else N=$3; fi

	if [ $1 = 'run' ] || [ $step_num = 7 ]; then
		if [ -e normalize ]; then rm -rf normalize; fi
		mkdir normalize
		$BIN/normalizeRTfile.pl -i $D -o ./normalize/D.norm.rt -m mean:vigintile2 -d 32 -l 32 -f 100
		$BIN/normalizeRTfile.pl -i $N -o ./normalize/N.norm.rt -m mean:vigintile2 -d 32 -l 32 -f 100
	fi

	##	8. calcErichment
	if [ $1 = 'run' ] || [ $step_num = 8 ]; then
        if [ -z $2 ]; then D=./normalize/D.norm.rt; else D=$2; fi
        if [ -z $3 ]; then N=./normalize/N.norm.rt; else N=$3; fi

		if [ -e enrichment ]; then rm -rf enrichment; fi
		mkdir enrichment
		$BIN/calcEnrich.pl -f $N -b $D -o ./enrichment/rRNA.erichment -w factor5:scaling1 -y 10 -x 0.25 -e complex
	fi

	##	9. Filter Erichment
	if [ $1 = 'run' ] || [ $step_num = 9 ]; then
		if [ -z $2 ]; then ERICH=./enrichment/rRNA.erichment; else ERICH=$2; fi
		if [ -z $3 ]; then OUT=rRNA.icSHAPE; else OUT=$3; fi

		$BIN/filterEnrich.pl -i $ERICH -o rRNA.icSHAPE -t 100 -T 2 -s 5 -e 30
	fi
fi
###########################################

###########################################                                  
##	icSAHPE Value to structure values
if [ $1 = 'cov' ]; then
	#cat rRNA.icSHAPE | sed 's/\t/\n/g' | sed -e '1,3d' | sed 's/NULL/0.000/g' > rRNA.varna
	if [ -z $2 ]; then FILE=rRNA.icSHAPE; else FILE=$2; fi
	if [ -z $3 ]; then FA=./DATA/fa/18SrRNA.fa; else FA=$3; fi
	if [ -z $4 ]; then DOTB=./DATA/fa/18SrRNA.dotbrak; else DOTB=$4; fi
	cat $FILE | awk '
	{split($0, Arr, "\t"); 
		i = 4
		while(i <= length(Arr))
		{
			if(Arr[i] == "NULL")
				print "0.000"
			else
				print Arr[i]
			i += 1
		}
	}
	' > rRNA.icSHAPE.varna
	cat $FA | grep -v "^>" | sed 's/\n//g'
	cat $DOTB | grep -v "^>" | sed 's/\n//g'
fi
###########################################
