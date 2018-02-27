#!/bin/bash

## Set parameters                                                                                                                                 
step=$1   #[PREP/INDIV/GROUP/GENO/VCF/INFO]
te=$2
DIR=$3
#DIR="/proj/sens2017106/nobackup/diana/melt/runsplit/" 
if [ ! $step ] || [ ! $te ] || [ ! -e $DIR ]; then
    echo "Usage: bash runMeltSplitSbatch.sh <STEP> <TE> <workDIR>"
    echo "STEP=PREP/INDIV/GROUP/GENO/VCF/INFO"
    echo "TE=Transposable Element, often: ALU/LINE1/SVA/(HERVK)"
    echo "workDIR: full path to where melt split has been setup by perl createMELTsplitSetup.pl"
    exit
fi

output=$DIR"intermediate/"
outIndiv=$output$te"DISCOVERY"
outGroup=$output$te"GROUP"
outGeno=$output$te"GENO"
results=$DIR"results/"

if [ $step = "PREP" ]; then
    for file in $DIR"01_PreProcess/ALL/sbatch/"*.sbatch
    do
        filebase=$(basename "$file")
        sample=${filebase%.*}
	dictfile=$DIR"bam/"$sample".dict"
	if [ -s $dictfile ]; then
	    echo "Sample $sample has been preprocessed. Remove $sample.dict to rerun"
	else
	    echo sbatch $file
	fi
    done

elif [ $step = "INDIV" ]; then
    adir=$DIR"02_IndivAnalysis/$te/sbatch"
    if [ ! -e $adir ]; then
	echo "Transposable element $te has not been setup, Run: perl createMELTsplitSetup.pl $te"
	exit
    fi
    for file in $adir/*.sbatch
    do
	filebase=$(basename "$file")
	sample=${filebase%.*}
	prepfile=$DIR"bam/"$sample".bam"
	if [ ! -s $prepfile ]; then
	    echo "File needs to be PreProcessed first ( $sample )"
	else
	    if [ -s $outIndiv"/"$sample"."$te".tmp.bed" ]; then
		echo "IndivAnalysis seems to have been run for $file"
	    else
		echo sbatch $file
	    fi
	fi
    done

elif [ $step = "GROUP" ]; then
    adir=$DIR"02_IndivAnalysis/$te/sbatch"
    for file in $adir/*.sbatch
    do
        filebase=$(basename "$file")
        sample=${filebase%.*}
	ires=$outIndiv"/"$sample"."$te".tmp.bed"
        if [ ! -s $ires ]; then
	    echo "Sample $sample has not completed IndivAnalysis"
	    exit
	fi
    done
    groupJob=$DIR"03_GroupAnalysis/"$te"/sbatch/GroupAnalysis.sbatch"
    echo sbatch $groupJob

elif [ $step = "GENO" ]; then
    resultsGroup=$outGroup"/"$te".pre_geno.tsv"
    if [ ! -s $resultsGroup ]; then
	echo "Need to run GroupAnalysis first"
    else
	for file in $DIR"04_Genotype/$te/sbatch/"*.sbatch
	do
	    filebase=$(basename "$file")
            sample=${filebase%.*}
            resfile=$outGeno"/"$sample.$te".tsv"
            indivfile=$outIndiv"/"$sample.$te".tmp.bed"
	    if [ -s $resfile ]; then
		echo "$resfile exists, skip (or rerun?)"
            elif [ ! -s $indivfile ]; then
		echo "Sample $sample has not been included in Indiv"
            else
		echo sbatch $file
	    fi
	done
    fi

elif [ $step = "VCF" ]; then
    adir=$DIR"04_Genotype/$te/sbatch"
    for file in $adir/*.sbatch
    do
        filebase=$(basename "$file")
        sample=${filebase%.*}
        gres=$outGeno"/"$sample"."$te".tsv"
        if [ ! -s $gres ]; then
            echo "Sample $sample has not completed Genotype"
            exit
        fi
    done
    groupVcf=$DIR"05_MakeVcf/"$te"/sbatch/MakeVcf.sbatch"
    echo sbatch $groupVcf

elif [ $step = "INFO" ]; then
    echo Preprocess: java -Xmx6G -jar $MELT_HOME"/MELT.jar" Preprocess bamfile ref
    echo Indiv java -Xmx6G -jar $MELT_HOME"/MELT.jar" IndivAnalysis -h ref -l bamfile -t te_zip -w $outIndiv -c covX -r readlength
    echo Group java -Xmx6G -jar $MELT_HOME"/MELT.jar" GroupAnalysis -h ref -l $outIndiv -t te_zip -w $outGroup -n genes_bed -r readlength
    echo Genotype java -Xmx6G -jar $MELT_HOME"/MELT.jar" Genotype -h ref -l bamfile -t te_zip -w $outGeno -p $outGroup -e insertsize
    echo MakeVCF java -Xmx6G -jar $MELT_HOME"/MELT.jar" MakeVCF -h ref -f tsv_list -t te_zip -w $outGeno -p $outGroup -o $results
else
    echo "Need to specify which step to run: INDIV/GROUP/GENO/VCF/INFO"
fi
