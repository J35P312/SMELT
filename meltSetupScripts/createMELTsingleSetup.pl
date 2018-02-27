#!/usr/bin/perl -w
use strict;

##################################################################################################################### 
# Name: createMELTsingleSetup.pl
# Desc: creates directories and sbatch scripts for each job that has to be run for MELTsplit for a given set of bam files
# Usage: perl createMELTsingleSetup.pl 

#my $meilist = $ARGV[0];
#($meilist || die "Usage perl createMELTsingleSetup.pl <MEI_loc_file.txt>\n");

#### USER DEFINED VARIABLES ##########################################################################################

# MELT settings
my $meltver = "2.0.2/bianca";
my $MELT_HOME="/sw/apps/bioinfo/MELT/".$meltver;
our $meltjar="java -Xmx6G -jar $MELT_HOME/MELT.jar";
our $genes_bed=$MELT_HOME."/add_bed_files/1KGP_Hg19/hg19.genes.bed";

#Uppmax project
our $project="sens2017106";
# Reference (has to be same as the one used to create bam files)
our $ref="/proj/".$project."/nobackup/private/wabi/melt/data/ref/hs37d5.fa";
our $excludeChrs="hs37d5/NC_007605"; #exclude decoy chrs, if any
# Directory where directories should be created
my $wdir = "/proj/".$project."/nobackup/diana/melt/runsingle/";
# File with paths to TE_zip files and optionally priors
my $meilist = $wdir . "mei_list_priors1kg.txt";
# File with full path to bam files to include in analysis
my $sample_file=$wdir."00_samples.txt";
# File with sample-id coverage insersize (this is optional, if values are similar for all bam files defaults can be used)
my $sample_stats = $wdir . "00_sample_stats.txt";
our $covX=30;
our $readlength=150;
our $insertsize=350;

######################################################################################################################

# Directories created by script in $wdir
our $bamlinkdir = $wdir."bam/";
my $outdir = $wdir."intermediate";
my $results = $wdir."results";

#### RUN MELT single (change times if needed) #########################################################################

unless(&checkInput($ref,$wdir,$meilist,$sample_file)) {exit;}
print "Set up for MELT single\nWorking dir: $wdir\nSamples:$sample_file\nMEI-models:$meilist\n";

&createDirStruct($outdir,$results,$bamlinkdir);
our @samples=&getSamples($sample_file);
our %sample_stats = &getSampleStats($sample_stats);
my @bamfiles = &linkBams();
#&preProcess("5:00:00");
&singleMELT("4:00:00");

#### SUBROUTINES #####################################################################################################

sub singleMELT {
    my ($time) = @_;
    my $pdir = "02_Single";
    my $dir_ok=&createDirectories($wdir.$pdir);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
        foreach my $bamfile (@bamfiles) {
            chomp $bamfile;
	    unless($bamfile =~ /.bam$/) {next;}
            my @path=split(/\//,$bamfile);
            my ($sample) = ($bamfile =~ /(.+)\.bam/);
            my $name="MELT_Single";
            my $cov = $covX; my $x = 0;
            #Use calculated value for covX if it exists    
            if(defined($sample_stats{$sample})) {
                ($cov,$x) = split(/\s/,$sample_stats{$sample});
            }
	    my $outSample = $outdir."/out_".$sample;
            my $cmd = "$meltjar Single -h $ref -b $excludeChrs -l $bamlinkdir$bamfile -n $genes_bed -t $meilist -w $outSample -c $covX -r $readlength -e $insertsize";
	    my @cmds=();
            push(@cmds,$cmd);
            &createSbatch($pdir,$sample,$time,$name,@cmds);
        }
    }
}


sub preProcess {
    my ($time) = @_;
    my $pdir = "01_PreProcess";
    my $te2 = "ALL";
    my $dir_ok=&createDirectories($pdir,$te2);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
	foreach my $bamfile (@bamfiles) {
	    chomp $bamfile;
	    unless($bamfile =~ /.bam$/) {next;}
	    my ($sample) = ($bamfile =~ /(.+)\.bam/);
	    my $name="MELT_PreP";
	    my $cmd = "$meltjar Preprocess $bamlinkdir$bamfile $ref";
	    my @cmds=();
	    push(@cmds,$cmd);
	    &createSbatch($pdir,$sample,$te2,$time,$name,@cmds);
	}
    }
}

sub linkBams {
    my @linkedbamfiles=();
    foreach my $bamfile (@samples) {
        chomp $bamfile;
	unless($bamfile =~ /.bam$/ && -s $bamfile) {next;}
        my @path=split(/\//,$bamfile);
        my $filebase = pop(@path);
        my ($sample) = ($filebase =~ /(.+)\.bam/);
        unless($sample =~ /[A-Z]/i) {print "Not correct format $bamfile\n"; next;}
        unless(-f $bamlinkdir.$filebase) {`ln -s $bamfile $bamlinkdir`;}
        #get bam index as well, either from src-dir or by samtools index, must be filebase.bam.bai                                                                
        push(@linkedbamfiles,$filebase);
        my $out_index = $bamlinkdir.$filebase.".bai";
        if(! -f $out_index) {
            my $index = $bamfile;
            $index =~ s/\.bam$/\.bai/;
            if(-f $index) {`ln -s $index $out_index`;}
        }
    }
    return(@linkedbamfiles);
}


sub createSbatch {
    my ($pdir,$sample,$time,$name,@commands) = @_;
    my $sbatch_file = $wdir.$pdir."/sbatch/".$sample.".sbatch";
    my $std_err_file = $wdir.$pdir."/std_err/".$sample.".err";
    my $std_out_file = $wdir.$pdir."/std_out/".$sample.".out";
    open(OUT,">$sbatch_file") || die "cannot open $sbatch_file";
    print OUT "#!/bin/bash\n";
    print OUT "#SBATCH -A ".$project."\n";
    print OUT "#SBATCH -p core\n";
    print OUT "#SBATCH -n 1\n";
    print OUT "#SBATCH -t $time\n";
    print OUT "#SBATCH -J ".$name."\n";
    print OUT "#SBATCH -o ".$std_out_file."\n";
    print OUT "#SBATCH -e ".$std_err_file."\n\n";
    print OUT "module load bioinfo-tools\nmodule load MELT/2.0.2\n";
    print OUT "module load bowtie2/2.3.3.1\n\n";
    foreach my $cmd (@commands) {
	chomp $cmd;
	print OUT "$cmd\n";
    }
    close OUT;
}

sub getSamples {
    my ($sfile) = @_;
    open(SIN,$sfile) || die "No sample file $sfile\n";
    my @sample_list=<SIN>;
    close SIN;
    return(@sample_list);
}

sub getSampleStats {
    my ($sfile) = @_;
    my %stats=();
    unless(-f $sfile) {return(%stats);}

    open(SIN,$sfile) || die "No sample file $sfile\n";
    while(<SIN>) {
	chomp;
	if(/^Sample/) {next;} #header
	my @line = split(/\s/);
	$stats{$line[0]} = "$line[1] $line[2]";
    }
    close SIN;
    return(%stats);
}

sub createDirectories {
    my ($dir_name) = @_;
    if(-e $dir_name) {
	return("Dir $dir_name exists, assuming has been created successfully");
    }
    else {
	`mkdir $dir_name`;
	`mkdir $dir_name"/std_err"`;
	`mkdir $dir_name"/std_out"`;
	`mkdir $dir_name"/sbatch"`;
	return("Successfully created $dir_name");
    }
}

sub createDirStruct {
    my ($tmpout,$resout,$bamout) = @_;
    unless(-e $tmpout) {
	`mkdir $tmpout`;
    }
    unless(-e $resout) {
	`mkdir $resout`;
	`mkdir $resout/vcf`;
    }
    unless(-e $bamout) {
        `mkdir $bamout`;
    }
}

sub checkInput {
    my ($refloc,$workdir,$meiloc,$samples) = @_;
    unless(-f $refloc) {
        print "No reference ( $refloc )\n";
        return(0);
    }
    unless(-e $workdir) {
	print "Dir $wdir does not exist\n";
        #`mkdir $workdir`;                                                                                                                                        
	return(0);
    }
    unless(-s $meiloc) {
        print "No MEI location file $meiloc\n";
        return(0);
    }
    unless(-s $samples) {
	print "No samples or sample file $samples\n";
        return(0);
    }
    return(1);
}


#"$meltjar Preprocess $bamfile $ref"                                                                                           
#"$meltjar IndivAnalysis -h $ref -l bamfile -t $te_zip -w $outIndiv -c $covX -r $readlength"                                                                        
#"$meltjar GroupAnalysis -h $ref -l $outIndiv -t $te_zip -w $outGroup -n $genes_bed -r $readlength"                                                                 
#"$meltjar Genotype -h $ref -l bamfile -t $te_zip -w $outGeno -p $outGroup -e $insertsize"                                                                          
#"$meltjar MakeVCF -h $ref -f $tsv_list -t $te_zip -w $outGeno -p $outGroup -o $results"              
