#!/bin/bash

############################################################################################################################################################################################################
###   	 ~BEFORE RUNNING~                                                                                                                                                                                ###
### Before running this script, please create a folder in $directory with $experiment_name with subdirectories '/raw_files/fastq' E.g. '~/Documents/Minion/My_Experiment/raw_files/fastq'. Then place    ###
### all generated fastq files in this folder. The same can be done for a fast5 folder in '/My_Experiment/raw_files/fast5'                                                                                ###
###                                                                                                                                                                                                      ###
############################################################################################################################################################################################################                                                                                                                                                                                                      

########### 	~Experiment details~     ########

experiment_name=
directory=~/Documents/Minion/$experiment_name

#################################################

### To create a folder structure for the current experiment with sample names. Provide a file in $directory with samples in barcoded order on separate lines called 'sample_names.txt'                   
### Sorting barcodes into folders will move all barcodes/trimmed barcodes into the created folder structure below. This should be done initially after demultiplexing and trimming, but should be turned 
### off when re-assembling.                                                                                                                                                                              

create_folder_structure=true
sort_barcodes_into_folders=true

#### To perform any of the following programs/assemblies, set variable to TRUE. Otherwise, set to FALSE.

demultiplexing=false
read_trimming=false
read_trim_length=2000
miniasm_assembly=true
indel_checking=false
auto_shutdown=false


if $demultiplexing; then
	cat $raw_minion_fastq/*fastq > $raw_minion_fastq/$experiment_name.cat.fastq
	porechop -i $raw_minion_fastq/$experiment_name.cat.fastq --discard_middle -t $threads -b $directory/barcodes
	rm $raw_minion_fastq/$experiment_name.cat.fastq
fi

if $read_trimming ; then
for sample in $directory/barcodes/BC*; do
	f=$(basename $sample)
	~/bin/Filtlong/bin/filtlong --min_length $read_trim_length --keep_percent 90 $directory/barcodes/$f > $directory/barcodes/trimmed.$f
	trimmed=true
done
fi

if $sort_barcodes_into_folders; then
	if $trimmed; then
		COUNTER=1
		for item in $directory/barcodes/trimmed.BC*; do
			f=$(basename $item)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			cp $directory/barcodes/$f $directory/samples/$location/reads/
			COUNTER=$[$COUNTER +1]
		done
	else
		COUNTER=1
		for item in $directory/barcodes/BC*; do
			f=$(basename $item)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			cp $directory/barcodes/$f $directory/samples/$location/reads/
			COUNTER=$[$COUNTER +1]
		done
	fi
fi
if $miniasm_assembly; then
	if $trimmed; then
		COUNTER=01
		for trimmed_sample in $directory/barcodes/trimmed.BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
			f=$(basename $trimmed_sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			mkdir $directory/samples/$location/miniasm_longread/

			~/bin/minimap2/minimap2 -t $threads -x ava-ont $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq | gzip -1 > $directory/samples/$location/miniasm_longread/minimap.gz

			~/bin/miniasm/miniasm -f $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.gz > $directory/samples/$location/miniasm_longread/miniasm.gfa

			awk '/^S/{print ">"$2"\n"$3}' $directory/samples/$location/miniasm_longread/miniasm.gfa > $directory/samples/$location/miniasm_longread/miniasm.fasta
			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/miniasm.fasta $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap.racon.paf

			racon -t $threads $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.racon.paf $directory/samples/$location/miniasm_longread/miniasm.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly.fasta

			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/consensus_assembly.fasta $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap_2.racon.paf

			racon -t $threads $directory/samples/$location/reads/trimmed.BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap_2.racon.paf $directory/samples/$location/miniasm_longread/consensus_assembly.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta

			assembly-stats $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta > assembly_stats_${location}.txt
			cp $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta $directory/samples/$location/miniasm_longread/${location}_miniasm_trimmed.fasta
		
			COUNTER=$[$COUNTER +1]
		done
	else
		COUNTER=01
		for sample in $directory/barcodes/BC*; do
			# $f not needed here, just generating it incase I want it later. This just uses number of barcodes to determine how many assemblies to produce.
			f=$(basename $sample)
			location=$(sed -n ${COUNTER}p $directory/sample_names.txt)
			mkdir $directory/samples/$location/miniasm_longread/

			~/bin/minimap2/minimap2 -t $threads -x ava-ont $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq | gzip -1 > $directory/samples/$location/miniasm_longread/minimap.gz

			~/bin/miniasm/miniasm -f $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.gz > $directory/samples/$location/miniasm_longread/miniasm.gfa

			awk '/^S/{print ">"$2"\n"$3}' $directory/samples/$location/miniasm_longread/miniasm.gfa > $directory/samples/$location/miniasm_longread/miniasm.fasta
			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/miniasm.fasta $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap.racon.paf

			racon -t $threads $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap.racon.paf $directory/samples/$location/miniasm_longread/miniasm.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly.fasta

			~/bin/minimap2/minimap2 -t $threads $directory/samples/$location/miniasm_longread/consensus_assembly.fasta $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq > $directory/samples/$location/miniasm_longread/minimap_2.racon.paf

			racon -t $threads $directory/samples/$location/reads/BC$(printf %02d $COUNTER).fastq $directory/samples/$location/miniasm_longread/minimap_2.racon.paf $directory/samples/$location/miniasm_longread/consensus_assembly.fasta > $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta
		
			assembly-stats $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta > assembly_stats_${location}.txt
			cp $directory/samples/$location/miniasm_longread/consensus_assembly2.fasta $directory/samples/$location/miniasm_longread/${location}_miniasm.fasta

			COUNTER=$[$COUNTER +1]
		done
	fi
fi
