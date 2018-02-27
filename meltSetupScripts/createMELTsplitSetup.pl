#!/usr/bin/perl -w
use strict;

##################################################################################################################### 
# Name: createMELTsplitSetup.pl
# Desc: creates directories and sbatch scripts for each job that has to be run for MELTsplit for a given set of bam files
# Usage: perl createMELTsplitSetup.pl <TE> where TE=ALU, LINE1 or SVA (HERVK)

my $te = $ARGV[0];
($te || die "Usage perl createPipelineStruct.pl <TE>\n");

#### USER DEFINED VARIABLES ##########################################################################################

# MELT settings
my $meltver = "2.0.2/bianca";
my $MELT_HOME="/sw/apps/bioinfo/MELT/".$meltver;
our $meltjar="java -Xmx6G -jar $MELT_HOME/MELT.jar";
our $te_zip=$MELT_HOME."/me_refs/1KGP_Hg19/".$te."_MELT.zip";
our $genes_bed=$MELT_HOME."/add_bed_files/1KGP_Hg19/hg19.genes.bed";

#Uppmax project
our $project="sens2017106";
# Reference (has to be same as the one used to create bam files)
our $ref="/proj/sens2017106/reference_material/fasta/human_g1k_v37.fasta";
our $excludeChrs="hs37d5/NC_007605"; #exclude decoy chrs, if any 
# Directory where directories should be created
my $wdir = "/proj/sens2017106/nobackup/diana/melt/swegen_analysis/output/";
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
our $outIndiv = $outdir . "/".$te."DISCOVERY";
our $outGroup = $outdir . "/".$te."GROUP";
our $outGeno = $outdir . "/".$te."GENO";
our $tsv_list=$outdir."/".$te."_tsv_list.txt";

#### RUN MELT split (change times if needed) #########################################################################

unless(&checkInput($ref,$wdir,$te_zip,$sample_file)) {exit;}

&createDirStruct($outdir,$results,$bamlinkdir);
our @samples=&getSamples($sample_file);
our %sample_stats = &getSampleStats($sample_stats);
my @bamfiles = &linkBams();
&preProcess("5:00:00");
&indivAnalysis("20:00:00");
&groupAnalysis("2:00:00");
&genotype("10:00:00");
&makeVcf("30:00");

#### SUBROUTINES #####################################################################################################

sub preProcess {
    my ($time) = @_;
    my $pdir = "01_PreProcess";
    my $te2 = "ALL";
    my $dir_ok=&createDirectories($wdir.$pdir,$te2);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
	foreach my $bamfile (@bamfiles) {
	    chomp $bamfile;
	    unless($bamfile =~ /.bam$/) {next;}
	    my ($sample) = ($bamfile =~ /(.+)\.bam/);
	    my $name="MELT_Prep";
	    my $cmd = "$meltjar Preprocess $bamlinkdir$bamfile $ref";
	    my @cmds=();
	    push(@cmds,$cmd);
	    &createSbatch($pdir,$sample,$te2,$time,$name,@cmds);
	}
    }
}

sub indivAnalysis {
    my ($time) = @_;
    my $pdir = "02_IndivAnalysis";
    my $dir_ok=&createDirectories($wdir.$pdir,$te);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
	foreach my $bamfile (@bamfiles) {
	    chomp $bamfile;
	    unless($bamfile =~ /.bam$/) {next;}
            my ($sample) = ($bamfile =~ /(.+)\.bam/);
	    my $name="MELT_IA";
	    my $cov = $covX; my $x = 0;
	    #Use calculated value for covX if it exists
	    if(defined($sample_stats{$sample})) {
		($cov,$x) = split(/\s/,$sample_stats{$sample});
	    }
	    my $cmd = "$meltjar IndivAnalysis -h $ref -b $excludeChrs -l $bamlinkdir$bamfile -t $te_zip -w \$TMPDIR/ -c $cov -r $readlength\ncp \$TMPDIR/* $outIndiv";
	    my @cmds=();
	    push(@cmds,$cmd);
	    &createSbatch($pdir,$sample,$te,$time,$name,@cmds);
	}
    }
}

sub groupAnalysis {
    my ($time) = @_;
    my $pdir = "03_GroupAnalysis";
    my $dir_ok=&createDirectories($wdir.$pdir,$te);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
	my $name="MELT_Gr";
	my $cmd = "$meltjar GroupAnalysis -h $ref -l $outIndiv -t $te_zip -w \$TMPDIR/ -n $genes_bed -r $readlength\ncp \$TMPDIR/* $outGroup";  
	my @cmds=();
	push(@cmds,$cmd);
	&createSbatch($pdir,"GroupAnalysis",$te,$time,$name,@cmds);
    }
}

sub genotype {
    my ($time) = @_;
    my $pdir = "04_Genotype";
    my $dir_ok=&createDirectories($wdir.$pdir,$te);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
	foreach my $bamfile (@bamfiles) {
	    chomp $bamfile;
	    unless($bamfile =~ /.bam$/) {next;}
            my ($sample) = ($bamfile =~ /(.+)\.bam/);
	    my $name="MELT_Ge";
	    my $ins = $insertsize; my $c = 0;
	    if(defined($sample_stats{$sample})) {
		($c,$ins) = split(/\s/,$sample_stats{$sample});
	    }
	    my $cmd = "$meltjar Genotype -h $ref -l $bamlinkdir$bamfile -t $te_zip -w \$TMPDIR/ -p $outGroup -e $ins\ncp \$TMPDIR/* $outGeno";
	    my @cmds=();
	    push(@cmds,$cmd);
	    &createSbatch($pdir,$sample,$te,$time,$name,@cmds);
	}
    }
}

sub makeVcf {
    my ($time) = @_;
    my $pdir = "05_MakeVCF";
    my $dir_ok=&createDirectories($wdir.$pdir,$te);
    print "$dir_ok\n";
    if($dir_ok =~ /^Successfully/) {
	my $name="MELT_Vcf";
	my $cmd1 = "printf".' "" > ' .$tsv_list."\n";
	$cmd1 .= "for tsvfile in $outGeno/*.".$te.".tsv\n";
	$cmd1 .= "do\n";
	$cmd1 .= "printf".' "%s\n" $tsvfile  >> ' .$tsv_list."\n";
	$cmd1 .= "done\n";
	
	my @cmds=();
	push(@cmds,$cmd1);
	my $cmd = "$meltjar MakeVCF -h $ref -f $tsv_list -t $te_zip -w $outGeno -p $outGroup -o $results"; 
	push(@cmds,$cmd);
	&createSbatch($pdir,"MakeVcf",$te,$time,$name,@cmds);
    }
}


sub createSbatch {
    my ($pdir,$sample,$tetmp,$time,$name,@commands) = @_;
    my $sbatch_file = $wdir.$pdir."/$tetmp/sbatch/".$sample.".sbatch";
    my $std_err_file = $wdir.$pdir."/".$tetmp."/std_err/".$sample.".err";
    my $std_out_file = $wdir.$pdir."/".$tetmp."/std_out/".$sample.".out";
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
        unless(-f $out_index) {
            my $index = $bamfile;
            $index =~ s/\.bam$/\.bai/;
            if(-f $index) {`ln -s $index $out_index`;}
        }
    }
    return(@linkedbamfiles);
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
    my ($cdir,$t_element) = @_;
    if(! -e $cdir) {
	`mkdir $cdir`; 
    }
    if(-e $cdir."/".$t_element) {
	return("Dir $cdir/$t_element exists, assuming has been created successfully");
    }
    else {
	`mkdir $cdir"/"$t_element`;
	`mkdir $cdir"/"$t_element"/std_err"`;
	`mkdir $cdir"/"$t_element"/std_out"`;
	`mkdir $cdir"/"$t_element"/sbatch"`;
	return("Successfully created $cdir/$t_element");
    }
}

sub createDirStruct {
    my ($tmpout,$resout,$bamout) = @_;
    unless(-e $tmpout) {
	`mkdir $tmpout`;
    }
    unless(-e $resout) {
	`mkdir $resout`;
	#`mkdir $resout/vcf`;
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
