# Author Ali Pirani and Dhatri Badri

configfile: "config/config.yaml"

import pandas as pd
import os
import gzip
import re

my_basedir = workflow.current_basedir
samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

#get name of ref genome
file_name = os.path.basename(config["reference_genome"])
REF_NAME = file_name.split('.')[0]
#print(REF_NAME)

#SHORTREADS = list(samples_df['sample_id'])
#print("Config:", config)

if not os.path.exists("results/" + PREFIX):
    try:
        os.makedirs("results/" + PREFIX)
    except OSError as e:
        print(f"Error creating directory: {e}")

# downsample reads
def downsample_reads(R1_file, R2_file, R1_out, R2_out, genome_size):
    
    R1_file = R1_file.pop()
    R2_file = R2_file.pop()
    R1_out = R1_out.pop()
    R2_out = R2_out.pop()

    gsize = genome_size.pop()

    print("Using Genome Size: %s to calculate coverage" % gsize)

    # Extract basic fastq reads stats with seqtk
    seqtk_check = "/nfs/esnitkin/bin_group/seqtk/seqtk fqchk -q3 %s > %s_fastqchk.txt" % (R1_file, R1_file)

    print(seqtk_check)
    
    try:
        os.system(seqtk_check)
    except sp.CalledProcessError: 
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % R1_file, 'rU') as file_open: 
        for line in file_open:
            if line.startswith('min_len'):
                line_split = line.split(';')
                min_len = line_split[0].split(': ')[1]
                max_len = line_split[1].split(': ')[1]
                avg_len = line_split[2].split(': ')[1]
            if line.startswith('ALL'):
                line_split = line.split('\t')
                total_bases = int(line_split[1]) * 2
    file_open.close()

    print('Average Read Length: %s' % avg_len)

    print('Total number of bases in fastq: %s' % total_bases)

    # Calculate original depth and check if it needs to be downsampled to a default coverage.
    ori_coverage_depth = int(total_bases / gsize)
    
    print('Original Covarage Depth: %s x' % ori_coverage_depth)

    #proc = sp.Popen(["nproc"], stdout=sp.PIPE, shell=True)
    #(nproc, err) = proc.communicate()
    #nproc = int(nproc.strip())

    if ori_coverage_depth > 100:
        # Downsample to 100
        factor = float(100 / float(ori_coverage_depth))
        #r1_sub = "/tmp/%s" % os.path.basename(R1_file)
        r1_sub = R1_out

        # Downsample using seqtk
        try:
            print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R1_file, factor, r1_sub)) 
            seqtk_downsample = "/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R1_file, factor, r1_sub) 
            os.system(seqtk_downsample)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if R2_file:
            #r2_sub = "/tmp/%s" % os.path.basename(R2_file)
            r2_sub = R2_out
            try:
                print("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R2_file, factor, r2_sub))  
                os.system("/nfs/esnitkin/bin_group/seqtk/seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R2_file, factor, r2_sub))  
            except sp.CalledProcessError:
                print('Error running seqtk for downsampling raw fastq reads.')
                sys.exit(1)
        else:
            r2_sub = "None"
    else:
        r1_sub = R1_file
        r2_sub = R2_file
        os.system("cp %s %s" % (R1_file, R1_out))
        os.system("cp %s %s" % (R2_file, R1_out))
        #if args.reverse_raw:
            #r2_sub = R2_file
        #else:
            #r2_sub = "None"
    #return r1_sub, r2_sub

def parse_bed_file(final_bed_unmapped_file):
    unmapped_positions_array = []
    with open(final_bed_unmapped_file, 'rU') as fp:
        for line in fp:
            line_array = line.split('\t')
            lower_index = int(line_array[1]) + 1
            upper_index = int(line_array[2]) + 1
            for positions in range(lower_index,upper_index):
                unmapped_positions_array.append(positions)
    only_unmapped_positions_file = final_bed_unmapped_file + "_positions"
    f1=open(only_unmapped_positions_file, 'w+')
    for i in unmapped_positions_array:
        p_string = str(i) + "\n"
        f1.write(p_string)
    return only_unmapped_positions_file

def remove_5_bp_snp_indel(raw_snp_vcf_file, raw_indel_vcf_file, output_file):
    remove_snps_5_bp_snp_indel_file_name = output_file
    #raw_snp_vcf_file + "_5bp_indel_removed.vcf"
    
    # Extract positions of indels
    indel_positions = set()
    with open(raw_indel_vcf_file, 'r') as indel_file:
        for line in indel_file:
            if not line.startswith('#'):
                line_array = line.split('\t')
                pos = int(line_array[1])
                indel_positions.add(pos)
    
    # Define range of positions to exclude
    exclude_positions = set()
    excluded_positions_list = []  # List to store excluded positions
    for indel_pos in indel_positions:
        for i in range(indel_pos - 5, indel_pos + 6):
            exclude_positions.add(i)
            
    # Print indel_positions and exclude_positions for sanity checks
    #print("Indel Positions:", indel_positions)
    #print("Exclude Positions:", exclude_positions)
    
    # Write filtered SNPs to new VCF file
    with open(remove_snps_5_bp_snp_indel_file_name, 'w') as filtered_file:
        with open(raw_snp_vcf_file, 'r') as snp_file:
            for line in snp_file:
                if line.startswith('#'):
                    filtered_file.write(line)
                else:
                    line_array = line.split('\t')
                    pos = int(line_array[1])
                    if pos in exclude_positions:
                        excluded_positions_list.append(pos)  # Store excluded positions
                    else:
                        filtered_file.write(line)
                    #if pos not in exclude_positions:
                        #filtered_file.write(line)
    
    # Print excluded positions for analysis
    #print("Excluded Positions:", excluded_positions_list)
    
    return remove_snps_5_bp_snp_indel_file_name

rule all:
    input:
        #trimmed_reads_forward=expand("results/{prefix}/{sample}/trimmomatic/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #trimmed_reads_reverse=expand("results/{prefix}/{sample}/trimmomatic/{sample}_R2_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #downsample_read_forward = expand("results/{prefix}/{sample}/downsample/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #downsample_read_reverse = expand("results/{prefix}/{sample}/downsample/{sample}_R2_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #aligned_reads = expand("results/{prefix}/{sample}/align_reads/{sample}_aln.sam", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #sorted_bam_reads= expand("results/{prefix}/{sample}/post_align/sorted_bam/{sample}_aln_sort.bam", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        dups_rmvd_sorted_bam_reads = expand("results/{prefix}/{sample}/post_align/sorted_bam_dups_removed/{sample}_final.bam", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #alignment_stats = expand("results/{prefix}/{sample}/stats/{sample}_alignment_stats.tsv", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #gatk_DoC = expand("results/{prefix}/{sample}/stats/{sample}_depth_of_coverage.sample_summary", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #unmapped_bam = expand("results/{prefix}/{sample}/bedtools/bedtools_unmapped/{sample}_unmapped.bed", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #bioawk_ref_size_file = expand("results/{prefix}/ref_genome_files/{ref_name}.size", prefix=PREFIX, ref_name=REF_NAME),
        #unmapped_bam_positions = expand("results/{prefix}/{sample}/bedtools/bedtools_unmapped/{sample}_unmapped.bed_positions", prefix=PREFIX, sample=SAMPLE, ref_name=REF_NAME),
        #bed_file = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME),
        #bedgraph_coverage = expand("results/{prefix}/{sample}/bedtools/bedgraph_coverage/{sample}.bedcov", prefix=PREFIX, sample=SAMPLE),
        #final_raw_gatk_vcf = expand("results/{prefix}/{sample}/gatk_varcall/{sample}_aln_mpileup_raw.vcf", prefix=PREFIX, sample=SAMPLE),
        #final_indel_vcf = expand("results/{prefix}/{sample}/gatk_varcall/{sample}_indel.vcf", prefix=PREFIX, sample=SAMPLE),
        #zipped_indel_vcf = expand("results/{prefix}/{sample}/gatk_varcall/{sample}_indel.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        final_raw_snp_vcf = expand("results/{prefix}/{sample}/samtools_varcall/{sample}_aln_mpileup_raw.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_final_raw_vcf = expand("results/{prefix}/{sample}/samtools_varcall/{sample}_aln_mpileup_raw.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        #filter_snp_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp.vcf", prefix=PREFIX, sample=SAMPLE),
        #filter_snp_final = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp_final.vcf", prefix=PREFIX, sample=SAMPLE),
        #filter_indel_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel.vcf", prefix=PREFIX, sample=SAMPLE),
        #filter_indel_final = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel_final.vcf", prefix=PREFIX, sample=SAMPLE),
        #zipped_filtered_snp_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp_final.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        #zipped_filtered_indel_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel_final.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        remove_snps_5_bp_snp_indel_file = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_5bp_indel_removed.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_5bp_indel_removed_snp_vcf_file = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_5bp_indel_removed.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        #freebayes_varcall = expand("results/{prefix}/{sample}/freebayes/{sample}_aln_freebayes_raw.vcf", prefix=PREFIX, sample=SAMPLE),
        #csv_summary_file = expand("results/{prefix}/{sample}/annotated_files/{sample}_ANN.csv", prefix=PREFIX, sample=SAMPLE),
        #annotated_vcf = expand("results/{prefix}/{sample}/annotated_files/{sample}_ANN.vcf", prefix=PREFIX, sample=SAMPLE),
        #zipped_freebayes_varcall = expand("results/{prefix}/{sample}/freebayes/{sample}_aln_freebayes_raw.vcf.gz", prefix=PREFIX, sample=SAMPLE)
        
        #unecesssary 
        #zipped_5bp_indel_removed_snp_vcf_file = expand("results/{prefix}/{sample}/remove_5_bp_snp_indel/{sample}_5bp_indel_removed.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        #zipped_filtered_snp_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp_final.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        #zipped_filtered_indel_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel_final.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        #zipped_indel_vcf = expand("results/{prefix}/{sample}/tabix/{sample}_filter_indel_final_zipped.gz", prefix=PREFIX, sample=SAMPLE),   
        #zipped_snp_vcf = expand("results/{prefix}/{sample}/tabix/{sample}_filtered_snp_final_zipped.gz", prefix=PREFIX, sample=SAMPLE),   
        #zipped_annotated_vcf = expand("results/{prefix}/{sample}/tabix/{sample}_ANN_zipped.gz", prefix=PREFIX, sample=SAMPLE)
        
# trims the raw fastq files to give trimmed fastq files
rule trim_reads:
    input:
        r1 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
    output:
        r1 = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz"),
        r2 = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz"), 
        r1_unpaired = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz"),
        r2_unpaired = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz")
    params:
        adapter_filepath=config["adaptor_filepath"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        threads = config["ncores"] # this uses bwa ncores
    log:
        trim_log = "logs/{prefix}/{sample}/trimmomatic/{sample}.log"
    #conda:
        #"envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log.trim_log}"

# downsamples trimmed fastq files
rule downsample:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz")
    output:
        outr1 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R1_trim_paired.fastq.gz",
        outr2 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R2_trim_paired.fastq.gz"
    params:
        gsize = config["genome_size"]
    log:
        downsample_log = "logs/{prefix}/{sample}/downsample/{sample}.log"
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})

# aligns trimmed fastq files using bwa to give sam file
rule align_reads:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/{wildcards.sample}_R2_trim_paired.fastq.gz")
    output:
        aligned_sam_out = temp(f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam")
        #touch
    params:
        #outdir = "results/{prefix}/{sample}/align_reads",
        num_cores = config["ncores"],
        ref_genome = config["reference_genome"],
        prefix = "{prefix}",
        sample = "{sample}",
        base_dir = my_basedir
        #prefix = "{sample}"
    log:
        bwa_log= "logs/{prefix}/{sample}/align_reads/{sample}.log"
    singularity:
        "docker://staphb/bwa:0.7.17"
    #conda:
        #"envs/bwa.yaml"
    shell:
        "./align_reads.sh {input.r1} {params.num_cores} {params.ref_genome} {input.r1} {input.r2} {output.aligned_sam_out}"
        #"""
        #split_field=$(python3 -c "from python_scripts.prepare_readgroup import prepare_readgroup; print(prepare_readgroup('{input.r1}'))") &&
        #bwa mem -M -R "$split_field" -t {params.num_cores} {params.ref_genome} {input.r1} {input.r2} > {output.aligned_sam_out}
        #"""
        #"""
        #mkdir -p {params.base_dir}/results/{params.prefix}/{params.sample}/align_reads/ &&
        #split_field=$(python3 -c "from python_scripts.prepare_readgroup import prepare_readgroup; print(prepare_readgroup('{params.base_dir}/{input.r1}'))") && 
        #bwa mem -M -R \"$split_field\" -t {params.num_cores} {params.ref_genome} {params.base_dir}/{input.r1} {params.base_dir}/{input.r2} > {params.base_dir}/{output.aligned_sam_out}
        #"""
        #"""
        #mkdir -p {params.base_dir}/results/{params.prefix}/{params.sample}/align_reads/ &&
        #split_field=$(zcat {input.r1} | awk 'NR==2{gsub("\t", "\\t"); print "@RG\\tID:"$1"\\tSM:"$1"\\tPL:ILLUMINA" }') &&
        #bwa mem -M -R "$split_field" -t {params.num_cores} {params.ref_genome} {input.r1} {input.r2} > {output.aligned_sam_out}
        #"""

# samclip and sort bam file
rule post_align_sam_to_bam:  
    input:
        aligned_sam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/align_reads/{wildcards.sample}_aln.sam")
    output:
        clipped_sam_out = temp(f"results/{{prefix}}/{{sample}}/post_align/samclip/{{sample}}_clipped.sam"),
        bam_out = temp(f"results/{{prefix}}/{{sample}}/post_align/aligned_bam/{{sample}}_aln.bam"),
        sorted_bam_out = temp(f"results/{{prefix}}/{{sample}}/post_align/sorted_bam/{{sample}}_aln_sort.bam")
    params:
        outdir_temp = "results/{prefix}/{sample}/post_align/sorted_bam/{sample}_aln_sort_temp",
        prefix = "{sample}",
        ref_genome= config["reference_genome"]
    wrapper:
        "file:python_scripts/sam_to_bam"

# remove duplicates and sort and index bam file with duplicates removed 
rule post_align_remove_duplicates:
    input:
        sorted_bam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam/{wildcards.sample}_aln_sort.bam")
    output:
        bam_duplicates_removed_out = temp(f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam"),
        dups_rmvd_sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sorted_bam_dups_removed/{{sample}}_final.bam"
    params:
        outdir_dups_removed = "results/{prefix}/{sample}/post_align/remove_duplicates",
        outdir = "results/{prefix}/{sample}/post_align/sorted_bam_dups_removed/",
        prefix = "{sample}"
    wrapper:
        "file:python_scripts/post_align_remove_duplicates"
    
# determine statistics of file
rule stats:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam")
    output:
        alignment_stats = f"results/{{prefix}}/{{sample}}/stats/{{sample}}_alignment_stats.tsv" 
    #conda:
        #"envs/samtools.yaml"
    singularity:
        "docker://staphb/samtools:1.19"
    shell:
        "samtools flagstat {input.index_sorted_dups_rmvd_bam} > {output.alignment_stats}" 

# determine coverage of bam file       
rule coverage_depth:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam")
    output:
        gatk_depthCoverage_summary = f"results/{{prefix}}/{{sample}}/stats/{{sample}}_depth_of_coverage.sample_summary"
    params:
        outdir = "results/{prefix}/{sample}/stats",
        ref_genome = config["reference_genome"], 
        prefix = "{sample}",
        intervals = config["bed_file"] # change to output from bedfile rule
    conda:
        "envs/gatk.yaml"
    shell:
        "gatk DepthOfCoverage -R {params.ref_genome} -O {params.outdir}/{params.prefix}_depth_of_coverage -I {input.index_sorted_dups_rmvd_bam} --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {params.intervals}"

# bedtools
rule bedtools:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam")
    output:
        unmapped_bed = f"results/{{prefix}}/{{sample}}/bedtools/bedtools_unmapped/{{sample}}_unmapped.bed"
    #conda:
        #"envs/bedtools.yaml"
    singularity:
        "docker://staphb/bedtools:2.31.1"
    shell:
        "bedtools genomecov -ibam {input.index_sorted_dups_rmvd_bam} -bga | awk '$4==0' > {output.unmapped_bed}"

# returns unmapped positions file    
rule parse_bed_file:
    input:
        unmapped_bed = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/bedtools/bedtools_unmapped/{wildcards.sample}_unmapped.bed")
    output:
        unmapped_bam_positions = f"results/{{prefix}}/{{sample}}/bedtools/bedtools_unmapped/{{sample}}_unmapped.bed_positions"
    run:
        parse_bed_file(input.unmapped_bed[0])

# prepared reference size and window files
rule prepare_references_files:
    #input:
        #reference_size_file = lambda wildcards: expand(f"results/{wildcards.prefix}/ref_genome_files/{wildcards.ref_name}.size")
    output:
        reference_size_file=f"results/{{prefix}}/ref_genome_files/{{ref_name}}.size",
        reference_window_file = f"results/{{prefix}}/ref_genome_files/{{ref_name}}.bed"
    params:
        ref_genome = config["reference_genome"]
    wrapper:
        "file:python_scripts/prepare_reference_files"

# created ref genome files 
#rule bioawk:
    #output:
        #reference_size_file=f"results/{{prefix}}/ref_genome_files/{{ref_name}}.size"
    #params:
        #ref_genome = config["reference_genome"]
    #conda:
        #"envs/bioawk.yaml"
    #shell:
        #"bioawk -c fastx '{{ print $name, length($seq) }}' < {params.ref_genome} > {output.reference_size_file}"

#create bed file from biowk
#rule create_bed_file:
    #input:
        #index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
        #reference_size_file = lambda wildcards: expand(f"results/{wildcards.prefix}/ref_genome_files/{wildcards.ref_name}.size")
    #output:
        #bedgraph_cov = f"results/{{prefix}}/{{sample}}/bedtools/bedgraph_coverage/{{sample}}.bedcov",
        #reference_window_file = f"results/{{prefix}}/ref_genome_files/{{ref_name}}.bed"
    #conda:
        #"envs/bedtools.yaml"
    #shell:
        #"bedtools makewindows -g {input.reference_size_file} -w 1000 > {output.reference_window_file}"

rule bedgraph_cov:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
        reference_window_file = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME)
    output:
        bedgraph_cov = f"results/{{prefix}}/{{sample}}/bedtools/bedgraph_coverage/{{sample}}.bedcov"
    singularity:
        "docker://staphb/bedtools:2.31.1"
    #conda:
        #"envs/bedtools.yaml"
    shell:
        "bedtools coverage -abam {input.index_sorted_dups_rmvd_bam} -b {input.reference_window_file} > {output.bedgraph_cov}"

# variant calling 
# gatk
# calling snp/indel and subset of variants using gatk
rule gatk_call_indel:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
    output:
        final_raw_vcf= f"results/{{prefix}}/{{sample}}/gatk_varcall/{{sample}}_aln_mpileup_raw.vcf",
        indel_file = f"results/{{prefix}}/{{sample}}/gatk_varcall/{{sample}}_indel.vcf",
        zipped_indel_vcf = f"results/{{prefix}}/{{sample}}/gatk_varcall/{{sample}}_indel.vcf.gz"
    params:
        #haplotype = config["haplotype_parameters"],
        ref_genome = config["reference_genome"]
    wrapper:
        "file:python_scripts/gatk_call_indel"

# calling snps with samtools
# variant_calling
rule bcftools_call_snps:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
    output:
        final_raw_vcf = f"results/{{prefix}}/{{sample}}/samtools_varcall/{{sample}}_aln_mpileup_raw.vcf",
        #zipped_final_raw_vcf = f"results/{{prefix}}/{{sample}}/samtools_varcall/{{sample}}_aln_mpileup_raw.vcf.gz"
        #final_raw_postalign_vcf = f"results/{{prefix}}/{{sample}}/samtools_varcall/{{sample}}_aln_mpileup_postalign_raw.vcf"
    params:
        ref_genome = config["reference_genome"], 
        #mpileup_params = config["mpileup_parameters"],
    singularity:
        "docker://staphb/bcftools:1.19"
    shell:
        """
        bcftools mpileup -f {params.ref_genome} {input.index_sorted_dups_rmvd_bam} | bcftools call -Ov -v -c -o {output.final_raw_vcf}
        """
    #wrapper:
        #"file:python_scripts/bcftools_call_snps"
        #module load Bioinformatics bcftools/1.12-g4b275e 

rule freebayes: 
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
    output:
        freebayes_varcall = f"results/{{prefix}}/{{sample}}/freebayes/{{sample}}_aln_freebayes_raw.vcf",
        zipped_freebayes_varcall = f"results/{{prefix}}/{{sample}}/freebayes/{{sample}}_aln_freebayes_raw.vcf.gz"
    params:
        ref_genome = config["reference_genome"], 
    conda:
        "envs/freebayes.yaml"
    #singularity:
        #"docker://staphb/freebayes:1.3.7",
        #"docker://staphb/htslib:1.19"
    shell:
        #"freebayes -f {params.ref_genome} {input.index_sorted_dups_rmvd_bam} > {output.freebayes_varcall}"
        """
        freebayes -f {params.ref_genome} {input.index_sorted_dups_rmvd_bam} > {output.freebayes_varcall}
        bgzip -c {output.freebayes_varcall} > {output.zipped_freebayes_varcall}
        tabix -p vcf -f {output.zipped_freebayes_varcall}
        """

rule hard_filter:
    input:
        final_raw_snp_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/samtools_varcall/{wildcards.sample}_aln_mpileup_raw.vcf"),
        final_raw_indel_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/gatk_varcall/{wildcards.sample}_indel.vcf")
    output:
        filter_snp_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp.vcf",
        filter_snp_final = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp_final.vcf",
        filter_indel_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel.vcf",
        filter_indel_final = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel_final.vcf",
        zipped_filtered_snp_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp_final.vcf.gz",
        zipped_filtered_indel_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel_final.vcf.gz"
    params:
        ref_genome = config["reference_genome"],
        dp_indel_filter = config["dp_indel_filter"],
        mq_indel_filter = config["mq_indel_filter"],
        qual_indel_filter = config["qual_indel_filter"],
        af_indel_filter = config["af_indel_filter"],
        dp_snp_filter = config["dp_snp_filter"],
        fq_snp_filter = config["fq_snp_filter"],
        mq_snp_filter = config["mq_snp_filter"],
        qual_snp_filter = config["qual_snp_filter"],
        af_snp_filter = config["af_snp_filter"]
    #conda:
        #"envs/gatk.yaml"
    wrapper:
        "file:python_scripts/hard_filter"

# remove snps with 5bp of an indel 
rule remove_5_bp_snp_indel:
    input:
        snp_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_snp_final.vcf"),
        indel_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_indel_final.vcf")
    output:
        remove_snps_5_bp_snp_indel_file = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_5bp_indel_removed.vcf",
        #zipped_indel_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_5bp_indel_removed.vcf.gz"
    run:
        remove_5_bp_snp_indel(input.snp_vcf[0], input.indel_vcf[0], output.remove_snps_5_bp_snp_indel_file)
    #wrapper:
        #"file:python_scripts/remove_5_bp_snp_indel"
    #shell:
        ##"""
        #remove_snps_5_bp_snp_indel_file=$(python3 -c "from python_scripts.remove_5_bp_snp_indel: import remove_5_bp_snp_indel; print(prepare_readgroup('{input.snp_vcf[0], input.indel_vcf[0], output.remove_snps_5_bp_snp_indel_file}'))")
        #"""


# create snpEff database from ref genome
rule install_annotate_snpEff:
    input:
        remove_snps_5_bp_snp_indel_file = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_5bp_indel_removed.vcf")
    output:
        #snpeff_db = f"data/{{ref_name}}/snpEffectPredictor.bin", ref_name=REF_NAME,
        csv_summary_file = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_ANN.csv",
        annotated_vcf = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_ANN.vcf"
    params:
        snpeff_parameters = config["snpeff_parameters"],
        #snpeff_path = config["snpEff_path"],
        #data_directory = config["data_dir"],
        #snpEff_config_file = config["snpEff_config_file"],
        snpEff_db = REF_NAME, # figure out how to call wildcards for ref_name 
        base_dir = my_basedir
    wrapper:
        "file:python_scripts/install_snpEff"

#rule annotate:
    #input:
        #remove_snps_5_bp_snp_indel_file = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/remove_5_bp_snp_indel/{wildcards.sample}_5bp_indel_removed.vcf"),
        #snpeff_db = expand("data/{ref_name}/snpEffectPredictor.bin", ref_name=REF_NAME)
    #output:
        #csv_summary_file = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_ANN.csv",
        #annotated_vcf = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_ANN.vcf"
    #params:
        #snpeff_parameters = config["snpeff_parameters"],
        #snpeff_path = config["snpEff_path"],
        #data_directory = config["data_dir"],
        #snpEff_db = lambda wildcards: REF_NAME,
        #snpEff_config_file = config["snpEff_config_file"]
    #shell:
        #"java -jar {params.snpeff_path} -csvStats {output.csv_summary_file} -dataDir {params.data_directory} {params.snpeff_parameters} -c {params.snpEff_config_file} {params.snpEff_db} {input.remove_snps_5_bp_snp_indel_file} > {output.annotated_vcf}" 

rule remove_5_bp_snp_indel_vcf:
    input:
        #indel_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_indel_final.vcf"),
        remove_snps_5_bp_snp_indel_file = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_5bp_indel_removed.vcf"),
        final_raw_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/samtools_varcall/{wildcards.sample}_aln_mpileup_raw.vcf"),
        #annotated_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/annotated_files/{wildcards.sample}_ANN.vcf")
    output:
        zipped_final_raw_vcf = f"results/{{prefix}}/{{sample}}/samtools_varcall/{{sample}}_aln_mpileup_raw.vcf.gz",
        zipped_remove_snps_5_bp_snp_indel_file = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_5bp_indel_removed.vcf.gz"
        #zipped_indel_vcf = f"results/{{prefix}}/{{sample}}/tabix/{{sample}}_filter_indel_final_zipped.gz",
        #zipped_snp_vcf = f"results/{{prefix}}/{{sample}}/tabix/{{sample}}_filtered_snp_final_zipped.gz", 
        #zipped_annotated_vcf = f"results/{{prefix}}/{{sample}}/tabix/{{sample}}_ANN_zipped.gz",
    #conda:
        #"envs/samtools.yaml"
    singularity:
        "docker://staphb/htslib:1.19"
    shell:
       """
       bgzip -c {input.remove_snps_5_bp_snp_indel_file} > {output.zipped_remove_snps_5_bp_snp_indel_file} 
       tabix -p vcf -f {output.zipped_remove_snps_5_bp_snp_indel_file}

       bgzip -c {input.final_raw_vcf} > {output.zipped_final_raw_vcf} && tabix -p vcf -f {output.zipped_final_raw_vcf}
       """