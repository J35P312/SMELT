//do not change these
params.input_dir="none"
params.help = false
params.step= false

//main script
if(params.help){
    println "SMELT: a Slurm wrapper for MELT"
    println "Usage: nextflow smelt.nf --step <STEP> --out <working_dir> --input <input_dir> --te <TE_zip_path> --ref {reference_fasta} -c smelt.conf"
    println "--step parameters"
    println "prep - run the MELT prep procedure on the bam files in the input_dir (step1)"
    println "indiv - run the indiv analysis on the extracted signals (step2)"
    println "group  - run the group step of the MELT split in input_dir (step3)"
    println "genotype  - run the MELT split analysis on the extracted signals in input_dir (step4)"
    println "makevcf  - run the MELT split analysis on the extracted signals in input_dir (step5)"
    println ""
    println "note! do not change the working_dir and input dir inbetween the analysis steps!"
    println ""
    println "other parameters:"
    println "--te - path to the transposable element zip file"

}else if (params.step == "prep"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    bam_files=Channel.fromPath("${params.folder}/*.bam").map{
        line ->
        ["${file(line).baseName}".replaceLast(/.bam/,""),file(line),line]
    }

    process prep{

        publishDir "${params.working_dir}/prep", mode: 'copy', overwrite: true
        scratch true
        errorStrategy 'ignore'


        cpus 1

        input:
            set ID, file(bam), bam_path from bam_files

        output:
            set "${bam.baseName}.bam.disc", "${bam.baseName}.bam.disc.bai", "${bam.baseName}.bam.fq" into prep_output 


        script:
        """
            samtools index ${bam}
            java -Xmx6G -jar ${params.melt} Preprocess ${bam} ${params.ref}
            mkdir ${params.working_dir}
            mkdir prep
            ln -s ${bam_path} ${params.working_dir}/${bam.baseName}.bam
        """

    }

}else if (params.step == "indiv"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    bam_files=Channel.fromPath("${params.folder}/*.bam").map{
        line ->
        ["${file(line).baseName}".replaceLast(/.bam/,""),file("${line}.disc"),file("${line}.disc.bai"),file("${line}.fq"),file(line)]
    }

    process indiv{

        publishDir "${params.working_dir}/indiv", mode: 'copy', overwrite: true
        scratch true
        errorStrategy 'ignore'

        cpus 1

        input:
            set ID, file(disc_bam), file(disc_index), file(fq),file(bam) from bam_files

        output:
            file "*.final.sorted.bam" into final_sorted
            file "*.final.sorted.bam.bai" into final_sorted_bai
            file "*.hum_breaks.sorted.bam" into hum_breaks
            file "*.hum_breaks.sorted.bam.bai" into hum_breaks_bai
            file "*.pulled.sorted.bam" into pulled_sorted
            file "*.pulled.sorted.bam.bai" into pulled_sorted_bai

        script:
        """
            mkdir tmp_indiv_out
            java -Xmx6G -jar ${params.melt} IndivAnalysis -l ${bam} -b {params.exclude} -h ${params.ref} -t ${params.te} -c ${params.cov} -r ${params.readlen} -w tmp_indiv_out
            cp tmp_indiv_out/* .
        """
    }

}else if (params.step == "group"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    process group{

        publishDir "${params.working_dir}/group", mode: 'copy', overwrite: true
        scratch true
        errorStrategy 'ignore'

        cpus 1
        
        input:
            file input_dir
        
        output:
            file "*" into outs

        script:


        """
        mkdir Group
        java -Xmx6G -jar ${params.melt} GroupAnalysis -h ${params.ref} -l ${input} -t ${params.te} -w Group -n ${params.genes} -r ${params.readlen}
        mv Group/* .
        rm -rf tmp

        """
    }

}else if (params.step == "genotype"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"

    bam_files=Channel.fromPath("${params.folder}/*.bam").map{
        line ->
        ["${file(line).baseName}".replaceLast(/.bam/,""),file(line)]
    }

    process genotype{

        publishDir "${params.working_dir}/genotype", mode: 'copy', overwrite: true
        scratch true
        errorStrategy 'ignore'


        cpus 1

        input:
            set ID, file(bam), bam_path from bam_files

        output:
            file "*" into genotype_output 


        script:
        """
            mkdir genotype_dir
            java -Xmx6G -jar ${params.melt} Genotype -h ${params.ref} -l ${bam} -t ${params.te} -w genotype_dir -p ${params.working_dir}/group/ -e ${insert_size}
        """

    }
}else if (params.step == "makevcf"){

    input_dir=file(params.input_dir)
    if(!input_dir.exists()) exit 1, "Error: The input directory was not found"



    process makevcf{

        publishDir "${params.working_dir}/genotype", mode: 'copy', overwrite: true
        scratch true
        errorStrategy 'ignore'


        cpus 1

        input:
            file input_dir

        output:
            file "*.vcf" into final_vcf 


        script:

        """
        printf "" > tsv_list.txt
        for tsvfile in ${params.working_dir}/genotype/*.tsv
            do
            printf "%s\n" $tsvfile  >> tsv_list.txt
        done

        mkdir vcf_output
        java -Xmx6G -jar ${params.melt} MakeVCF -h ${params.ref} -f tsv_list.txt -t ${params.te} -w ${params.working_dir}/genotype/ -p ${params.working_dir}/group -o vcf_output
        mv vcf_output/*vcf .
        """

    }

}else{
    println  "type: nextflow smelt.nf --help"
}
