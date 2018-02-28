# SMELT

SMELT is a Slurm wrapper around MELT. SMELT has been tested on the UPPMAX bianca cluster. With some luck, it might work on your cluster aswell!


# Install

Download and install MELT:

	http://melt.igs.umaryland.edu/manual.php

Edit the variables  in the createMELTsplitSetup.pl script, these are found at lines 15-44.

# Setup the output folder

	0: create a folder were the output files are stored, and update the $wdir parameter in the createMELTsplitSetup.pl script:
		mkdir path/working_dir

		and set : my $wdir = "path/working_dir"

	1: create a folder containing softlinks to each bam to analyse:
		
		mkdir path/working_dir/bam

	now create the softlinks:

		ln -s /your_bam_folder/*.ba* path/working_dir/bam

	2: create the 00_samples.txt file. This file contain the path to each bam file to analyse:

		ls  path/working_dir/bam/*bam > path/working_dir/00_samples.txt

	3: run the createMELTsplitSetup.pl. This script will create the other scripts and folders
	
		./createMELTsplitSetup.pl <TE>

	The script is run once for each <TE> type, such as ALU or LINE1


		./createMELTsplitSetup.pl ALU

		./createMELTsplitSetup.pl LINE1

		./createMELTsplitSetup.pl SVA

# Run the analysis

	Run the runMeltSplitSbatch.sh script. The script is  run once for each TE type and step:

		runMeltSplitSbatch.sh <STEP> <TE> <workDIR>
	
# Contributors

	Diana Ekman
	Jesper Eisfeldt
	Daniel Nilsson
