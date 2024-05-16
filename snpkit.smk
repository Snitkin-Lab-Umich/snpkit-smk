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
        os.system("cp %s %s" % (R2_file, R2_out))

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

def remove_5_bp_snp_indel(raw_snp_vcf_file, raw_indel_vcf_file, output_file, excluded_positions_file):
    remove_snps_5_bp_snp_indel_file_name = output_file
    
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
            excluded_positions_list.append(i)  # Store excluded positions
            
    # Print indel_positions and exclude_positions for sanity checks
    #print("Indel Positions:", indel_positions)
    #print("Exclude Positions:", exclude_positions)
    
    # Write excluded positions to the specified file
    with open(excluded_positions_file, 'w') as excluded_file:
        for pos in excluded_positions_list:
            excluded_file.write(str(pos) + '\n')

    # Write filtered SNPs to new VCF file
    with open(remove_snps_5_bp_snp_indel_file_name, 'w') as filtered_file:
        with open(raw_snp_vcf_file, 'r') as snp_file:
            for line in snp_file:
                if line.startswith('#'):
                    filtered_file.write(line)
                else:
                    line_array = line.split('\t')
                    pos = int(line_array[1])
                    #if pos in exclude_positions:
                        #excluded_positions_list.append(pos)  # Store excluded positions
                    #else:
                        #filtered_file.write(line)
                    if pos not in exclude_positions:
                        filtered_file.write(line)
    
    # Print excluded positions for analysis
    #print("Excluded Positions:", excluded_positions_list)
    
    return remove_snps_5_bp_snp_indel_file_name

rule all:
    input:
        trimmed_reads_forward = expand("results/{prefix}/{sample}/trimmomatic/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        trimmed_reads_reverse = expand("results/{prefix}/{sample}/trimmomatic/{sample}_R2_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        downsample_read_forward = expand("results/{prefix}/{sample}/downsample/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        downsample_read_reverse = expand("results/{prefix}/{sample}/downsample/{sample}_R2_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        aligned_reads = expand("results/{prefix}/{sample}/align_reads/{sample}_aln.sam", prefix=PREFIX, sample=SAMPLE),
        sorted_bam_reads= expand("results/{prefix}/{sample}/post_align/sorted_bam/{sample}_aln_sort.bam", prefix=PREFIX, sample=SAMPLE),
        dups_rmvd_sorted_bam_reads = expand("results/{prefix}/{sample}/post_align/sorted_bam_dups_removed/{sample}_final.bam", prefix=PREFIX, sample=SAMPLE),
        alignment_stats = expand("results/{prefix}/{sample}/stats/{sample}_alignment_stats.tsv", prefix=PREFIX, sample=SAMPLE),
        gatk_DoC = expand("results/{prefix}/{sample}/stats/{sample}_depth_of_coverage.sample_summary", prefix=PREFIX, sample=SAMPLE),
        unmapped_bam = expand("results/{prefix}/{sample}/bedtools/bedtools_unmapped/{sample}_unmapped.bed", prefix=PREFIX, sample=SAMPLE),
        bioawk_ref_size_file = expand("results/{prefix}/ref_genome_files/{ref_name}.size", prefix=PREFIX, ref_name=REF_NAME),
        unmapped_bam_positions = expand("results/{prefix}/{sample}/bedtools/bedtools_unmapped/{sample}_unmapped.bed_positions", prefix=PREFIX, sample=SAMPLE),
        bed_file = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME),
        bedgraph_coverage = expand("results/{prefix}/{sample}/bedtools/bedgraph_coverage/{sample}.bedcov", prefix=PREFIX, sample=SAMPLE),
        final_raw_gatk_vcf = expand("results/{prefix}/{sample}/gatk_varcall/{sample}_aln_mpileup_raw.vcf", prefix=PREFIX, sample=SAMPLE),
        final_indel_vcf = expand("results/{prefix}/{sample}/gatk_varcall/{sample}_indel.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_indel_vcf = expand("results/{prefix}/{sample}/gatk_varcall/{sample}_indel.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        final_raw_snp_vcf = expand("results/{prefix}/{sample}/samtools_varcall/{sample}_aln_mpileup_raw.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_final_raw_vcf = expand("results/{prefix}/{sample}/samtools_varcall/{sample}_aln_mpileup_raw.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        filter_snp_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp.vcf", prefix=PREFIX, sample=SAMPLE),
        filter_snp_final = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp_final.vcf", prefix=PREFIX, sample=SAMPLE),
        filter_indel_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel.vcf", prefix=PREFIX, sample=SAMPLE),
        filter_indel_final = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel_final.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_filtered_snp_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp_final.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        zipped_filtered_indel_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel_final.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        remove_snps_5_bp_snp_indel_file = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_5bp_indel_removed.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_5bp_indel_removed_snp_vcf_file = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_5bp_indel_removed.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        freebayes_varcall = expand("results/{prefix}/{sample}/freebayes/{sample}_aln_freebayes_raw.vcf", prefix=PREFIX, sample=SAMPLE),
        csv_summary_file = expand("results/{prefix}/{sample}/annotated_files/{sample}_ANN.csv", prefix=PREFIX, sample=SAMPLE),
        annotated_vcf = expand("results/{prefix}/{sample}/annotated_files/{sample}_ANN.vcf", prefix=PREFIX, sample=SAMPLE),
        zipped_freebayes_varcall = expand("results/{prefix}/{sample}/freebayes/{sample}_aln_freebayes_raw.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        zipped_filter_snp_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_snp.vcf.gz", prefix=PREFIX, sample=SAMPLE),
        zipped_filter_indel_vcf = expand("results/{prefix}/{sample}/filtered_vcf/{sample}_filter_indel.vcf.gz", prefix=PREFIX, sample=SAMPLE)
        
        
# trims the raw fastq files to give trimmed fastq files
rule trim_raw_reads:
    input:
        r1 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
    output:
        #r1 = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz"),
        #r2 = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz"), 
        #r1_unpaired = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz"),
        #r2_unpaired = temp(f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz")
        r1 = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz",
        r2 = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz", 
        r1_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz"
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
        threads = config["ncores"] 
    log:
        trim_log = "logs/{prefix}/{sample}/trimmomatic/{sample}.log"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log.trim_log}"

# downsamples trimmed fastq files
rule downsample_clean_reads:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz")
    output:
        outr1 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R1_trim_paired.fastq.gz",
        outr2 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R2_trim_paired.fastq.gz"
    params:
        gsize = config["genome_size"]
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize})

# aligns trimmed fastq files using bwa to give sam file
rule align_reads:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/downsample/{wildcards.sample}_R2_trim_paired.fastq.gz")
    output:
        #aligned_sam_out = temp(f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam")
        aligned_sam_out = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam"
    params:
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
    shell:
        "./align_reads.sh {input.r1} {params.num_cores} {params.ref_genome} {input.r1} {input.r2} {output.aligned_sam_out} &> {log.bwa_log}"

# samclip and sort bam file
rule removed_clipped_reads:  
    input:
        aligned_sam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/align_reads/{wildcards.sample}_aln.sam")
    output:
        #clipped_sam_out = temp(f"results/{{prefix}}/{{sample}}/post_align/samclip/{{sample}}_clipped.sam"),
        clipped_sam_out = f"results/{{prefix}}/{{sample}}/post_align/samclip/{{sample}}_clipped.sam",
        #bam_out = temp(f"results/{{prefix}}/{{sample}}/post_align/aligned_bam/{{sample}}_aln.bam"),
        bam_out = f"results/{{prefix}}/{{sample}}/post_align/aligned_bam/{{sample}}_aln.bam",
        #sorted_bam_out = temp(f"results/{{prefix}}/{{sample}}/post_align/sorted_bam/{{sample}}_aln_sort.bam")
        sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sorted_bam/{{sample}}_aln_sort.bam"
    params:
        outdir_temp = "results/{prefix}/{sample}/post_align/sorted_bam/{sample}_aln_sort_temp",
        prefix = "{sample}",
        ref_genome= config["reference_genome"]
    log:
        samclip_log = "logs/{prefix}/{sample}/post_align/{sample}_samclip.log",
        post_align_sam_to_bam_log= "logs/{prefix}/{sample}/post_align/{sample}.log"
    wrapper:
        "file:wrapper_functions/sam_to_bam"

# remove duplicates and sort and index bam file with duplicates removed 
rule post_align_remove_pcr_duplicates:
    input:
        sorted_bam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam/{wildcards.sample}_aln_sort.bam")
    output:
        #bam_duplicates_removed_out = temp(f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam"),
        bam_duplicates_removed_out = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam",
        dups_rmvd_sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sorted_bam_dups_removed/{{sample}}_final.bam"
    params:
        outdir_dups_removed = "results/{prefix}/{sample}/post_align/remove_duplicates",
        outdir = "results/{prefix}/{sample}/post_align/sorted_bam_dups_removed/",
        prefix = "{sample}"
    log:
        picard_dups_log = "logs/{prefix}/{sample}/post_align/{sample}_picard.log",
        post_align_remove_duplicates_log= "logs/{prefix}/{sample}/post_align/{sample}_samtools.log"
    wrapper:
        "file:wrapper_functions/post_align_remove_duplicates"
    
# determine statistics of file
rule alignment_stats:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam")
    output:
        alignment_stats = f"results/{{prefix}}/{{sample}}/stats/{{sample}}_alignment_stats.tsv" 
    singularity:
        "docker://staphb/samtools:1.19"
    log:
        "logs/{prefix}/{sample}/post_align/{sample}_stats.log"
    shell:
        "samtools flagstat {input.index_sorted_dups_rmvd_bam} > {output.alignment_stats} &> {log}" 

# determine coverage of bam file       
rule gatk_coverage_depth_statistics:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam")
    output:
        gatk_depthCoverage_summary = f"results/{{prefix}}/{{sample}}/stats/{{sample}}_depth_of_coverage.sample_summary"
    params:
        outdir = "results/{prefix}/{sample}/stats",
        ref_genome = config["reference_genome"], 
        prefix = "{sample}",
        intervals = config["bed_file"] 
    log:
        "logs/{prefix}/{sample}/stats/{sample}_coverage_depth.log"
    #conda:
        #"envs/gatk.yaml"
    singularity:
        "docker://broadinstitute/gatk:4.5.0.0"
    shell:
        "gatk DepthOfCoverage -R {params.ref_genome} -O {params.outdir}/{params.prefix}_depth_of_coverage -I {input.index_sorted_dups_rmvd_bam} --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {params.intervals} &> {log}"

# bedtools
rule bedtools_extract_coverage:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam")
    output:
        unmapped_bed = f"results/{{prefix}}/{{sample}}/bedtools/bedtools_unmapped/{{sample}}_unmapped.bed"
    log:
       "logs/{prefix}/{sample}/bedtools_unmapped/{sample}.log"
    singularity:
        "docker://staphb/bedtools:2.31.1"
    shell:
        "bedtools genomecov -ibam {input.index_sorted_dups_rmvd_bam} -bga | awk '$4==0' > {output.unmapped_bed}"

# returns unmapped positions file    
rule parse_bed_file_find_unmapped_regions:
    input:
        unmapped_bed = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/bedtools/bedtools_unmapped/{wildcards.sample}_unmapped.bed")
    output:
        unmapped_bam_positions = f"results/{{prefix}}/{{sample}}/bedtools/bedtools_unmapped/{{sample}}_unmapped.bed_positions"
    run:
        parse_bed_file(input.unmapped_bed[0])

# prepared reference size and window files
rule prepare_reference_windows:
    output:
        reference_size_file=f"results/{{prefix}}/ref_genome_files/{{ref_name}}.size",
        reference_window_file = f"results/{{prefix}}/ref_genome_files/{{ref_name}}.bed"
    params:
        ref_genome = config["reference_genome"]
    wrapper:
        "file:wrapper_functions/prepare_reference_files"

rule bedcoverage:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
        reference_window_file = expand("results/{prefix}/ref_genome_files/{ref_name}.bed", prefix=PREFIX, ref_name=REF_NAME)
    output:
        bedgraph_cov = f"results/{{prefix}}/{{sample}}/bedtools/bedgraph_coverage/{{sample}}.bedcov"
    singularity:
        "docker://staphb/bedtools:2.31.1"
    log:
        "logs/{prefix}/{sample}/bedgraph_coverage/{sample}_bedcov.log"
    shell:
        "bedtools coverage -abam {input.index_sorted_dups_rmvd_bam} -b {input.reference_window_file} > {output.bedgraph_cov} 2>&1 | tee {log}"
        
# variant calling 
# gatk
# calling snp/indel and subset of variants using gatk
rule gatk_call_indels:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
    output:
        final_raw_vcf= f"results/{{prefix}}/{{sample}}/gatk_varcall/{{sample}}_aln_mpileup_raw.vcf",
        indel_file = f"results/{{prefix}}/{{sample}}/gatk_varcall/{{sample}}_indel.vcf",
        zipped_indel_vcf = f"results/{{prefix}}/{{sample}}/gatk_varcall/{{sample}}_indel.vcf.gz"
    params:
        ref_genome = config["reference_genome"]
    log:
        "logs/{prefix}/{sample}/gatk/{sample}_gatk.log"
    wrapper:
        "file:wrapper_functions/gatk_call_indel"

# calling snps with samtools
# variant_calling
rule bcftools_call_snps:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
    output:
        final_raw_vcf = f"results/{{prefix}}/{{sample}}/samtools_varcall/{{sample}}_aln_mpileup_raw.vcf",
    log:
        "logs/{prefix}/{sample}/samtools_varcall/{sample}_bcftools_call_snps.log"
    params:
        ref_genome = config["reference_genome"]
    singularity:
        "docker://staphb/bcftools:1.19"
    shell:
        "bcftools mpileup -f {params.ref_genome} {input.index_sorted_dups_rmvd_bam} | bcftools call -Ov -v -c -o {output.final_raw_vcf} &> {log}"

rule freebayes_allele_frequency: 
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam_dups_removed/{wildcards.sample}_final.bam"),
    output:
        freebayes_varcall = f"results/{{prefix}}/{{sample}}/freebayes/{{sample}}_aln_freebayes_raw.vcf",
        #zipped_freebayes_varcall = f"results/{{prefix}}/{{sample}}/freebayes/{{sample}}_aln_freebayes_raw.vcf.gz"
    log:
        "logs/{prefix}/{sample}/freebayes/{sample}_freebayes.log"
    params:
        ref_genome = config["reference_genome"], 
    #conda:
    #    "envs/freebayes.yaml"
    singularity:
        "docker://staphb/freebayes:1.3.7"
    shell:
        """
        freebayes -f {params.ref_genome} {input.index_sorted_dups_rmvd_bam} > {output.freebayes_varcall} &> {log}
        """
        #bgzip -c {output.freebayes_varcall} > {output.zipped_freebayes_varcall} && tabix -p vcf -f {output.zipped_freebayes_varcall} &> {log}
        

rule variant_hard_filter:
    input:
        final_raw_snp_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/samtools_varcall/{wildcards.sample}_aln_mpileup_raw.vcf"),
        final_raw_indel_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/gatk_varcall/{wildcards.sample}_indel.vcf")
    output:
        filter_snp_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp.vcf",
        filter_snp_final = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp_final.vcf",
        filter_indel_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel.vcf",
        filter_indel_final = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel_final.vcf",
        zipped_filtered_snp_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp_final.vcf.gz",
        zipped_filtered_indel_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel_final.vcf.gz",
        zipped_filter_snp_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_snp.vcf.gz",
        zipped_filter_indel_vcf = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_filter_indel.vcf.gz",
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
    log:
        gatk_snp = "logs/{prefix}/{sample}/gatk/{sample}_filtered_snps.log",
        gatk_indels = "logs/{prefix}/{sample}/gatk/{sample}_filtered_indels.log"
    wrapper:
        "file:wrapper_functions/hard_filter"

# remove snps with 5bp of an indel 
rule remove_5_bp_snp_flanking_to_indels:
    input:
        snp_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_snp_final.vcf"),
        indel_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_indel_final.vcf")
    output:
        remove_snps_5_bp_snp_indel_file = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_5bp_indel_removed.vcf",
        excluded_positions_file = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_no_proximate_snp.vcf_positions_array"
    run:
        remove_5_bp_snp_indel(input.snp_vcf[0], input.indel_vcf[0], output.remove_snps_5_bp_snp_indel_file, output.excluded_positions_file)

# create snpEff database from ref genome
rule snpEff_annotation:
    input:
        remove_snps_5_bp_snp_indel_file = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_5bp_indel_removed.vcf"),
        filter_snp_final = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_snp_final.vcf"),
        filter_indel_final = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_filter_indel_final.vcf")
    output:
        csv_summary_file = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_ANN.csv",
        annotated_vcf = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_ANN.vcf",
        annotated_filter_snp_final = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_filter_snp_final_ANN.vcf",
        annotated_filter_indel_final = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_filter_indel_final_ANN.vcf",
        summary_file_filter_snp_final = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_filter_snp_final_ANN.csv",
        summary_file_filter_indel_final = f"results/{{prefix}}/{{sample}}/annotated_files/{{sample}}_filter_indel_final_ANN.csv"
    params:
        snpeff_parameters = config["snpeff_parameters"],
        snpEff_db = REF_NAME, # figure out how to call wildcards for ref_name 
        base_dir = my_basedir
    log:
        "logs/{prefix}/{sample}/snpEff/{sample}_snpEff.log"
    wrapper:
        "file:wrapper_functions/install_snpEff"

rule tabix_index_vcf:
    input:
        remove_snps_5_bp_snp_indel_file = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/filtered_vcf/{wildcards.sample}_5bp_indel_removed.vcf"),
        final_raw_vcf = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/samtools_varcall/{wildcards.sample}_aln_mpileup_raw.vcf"),
        freebayes_varcall = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/freebayes/{wildcards.sample}_aln_freebayes_raw.vcf",)
    output:
        zipped_final_raw_vcf = f"results/{{prefix}}/{{sample}}/samtools_varcall/{{sample}}_aln_mpileup_raw.vcf.gz",
        zipped_remove_snps_5_bp_snp_indel_file = f"results/{{prefix}}/{{sample}}/filtered_vcf/{{sample}}_5bp_indel_removed.vcf.gz",
        zipped_freebayes_varcall = f"results/{{prefix}}/{{sample}}/freebayes/{{sample}}_aln_freebayes_raw.vcf.gz"
    singularity:
        "docker://staphb/htslib:1.19"
    shell:
       """
       bgzip -c {input.remove_snps_5_bp_snp_indel_file} > {output.zipped_remove_snps_5_bp_snp_indel_file} 
       tabix -p vcf -f {output.zipped_remove_snps_5_bp_snp_indel_file}

       bgzip -c {input.final_raw_vcf} > {output.zipped_final_raw_vcf} && tabix -p vcf -f {output.zipped_final_raw_vcf}

       bgzip -c {input.freebayes_varcall} > {output.zipped_freebayes_varcall} && tabix -p vcf -f {output.zipped_freebayes_varcall} 
       """
