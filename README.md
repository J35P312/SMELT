# SMELT

SMELT is a nextflow wrapper around MELT.

# Install

Download and install MELT:

	http://melt.igs.umaryland.edu/manual.php

Then install nextflow:

	https://www.nextflow.io/

Edit the path variables and settings in the config file:

	nano smelt.conf

now you are ready to go!

#run melt split analysis

The script is run once for each step:

	./nextflow smelt.nf -c smelt.conf --step prep --working_dir <absolute/path/to/output/folder/> --input_dir <absolute/path/to/folder/containing/bam_files> --te <transposable_element_path>

	./nextflow smelt.nf -c smelt.conf --step indiv --working_dir <absolute/path/to/output/folder/> --input_dir <absolute/path/to/folder/containing/bam_files> --te <transposable_element_path>

	./nextflow smelt.nf -c smelt.conf --step group --working_dir <absolute/path/to/output/folder/> --input_dir <absolute/path/to/folder/containing/bam_files> --te <transposable_element_path>

	./nextflow smelt.nf -c smelt.conf --step genotype --working_dir <absolute/path/to/output/folder/> --input_dir <absolute/path/to/folder/containing/bam_files> --te <transposable_element_path>

	./nextflow smelt.nf -c smelt.conf --step makevcf --working_dir <absolute/path/to/output/folder/> --input_dir <absolute/path/to/folder/containing/bam_files> --te <transposable_element_path>

the output will be stored in the --working_dir folder, and the script will analyse all bam files in the --input_dir folder. The --te file is the mobile-element zip. These commands needs to be run separately for each mobile element type.
It is safest to select a separate working_dir for each mobile element type.

# Run melt single analysis
	
	./nextflow smelt.nf -c smelt.conf --step single --working_dir <absolute/path/to/output/folder/> --input_dir <absolute/path/to/folder/containing/bam_files> --te <transposable_element_path>

the output will be stored in the --working_dir folder, and the script will analyse all bam files in the --input_dir folder. The --te file is either the mobile-element zip file or a text file containing absolute path
to multiple such files (see the MELT documentation for more info) 

# Contributors

	Diana Ekman
	Jesper Eisfeldt
	Daniel Nilsson
