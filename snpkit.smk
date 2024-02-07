# Author Ali Pirani and Dhatri Badri

configfile: "config/config.yaml"

import pandas as pd
import os
import gzip
import re

samples_df = pd.read_csv(config["samples"])
SAMPLE = list(samples_df['sample_id'])

PREFIX = config["prefix"]

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

# Define a function to get all reference genome names
#def get_all_ref_names(config):
    #return [get_ref_name(ref_genome_path) for ref_genome_path in config["reference_genome"]]
# Get all reference genome names
#all_ref_names = get_all_ref_names(config)

def get_ref_name(ref_genome_path):
    file_name = os.path.basename(ref_genome_path)
    ref_name = file_name.split('.')[0]
    return ref_name

rule all:
    input:
        trimmed_reads_forward=expand("results/{prefix}/{sample}/trimmomatic/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        trimmed_reads_reverse=expand("results/{prefix}/{sample}/trimmomatic/{sample}_R2_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        downsample_read_forward = expand("results/{prefix}/{sample}/downsample/{sample}_R1_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        downsample_read_reverse = expand("results/{prefix}/{sample}/downsample/{sample}_R2_trim_paired.fastq.gz", prefix=PREFIX, sample=SAMPLE),
        aligned_reads = expand("results/{prefix}/{sample}/align_reads/{sample}_aln.sam", prefix=PREFIX, sample=SAMPLE),
        clipped_sam_files = expand("results/{prefix}/{sample}/post_align/samclip/{sample}_clipped.sam", prefix=PREFIX, sample=SAMPLE),
        aligned_bam_files = expand("results/{prefix}/{sample}/post_align/aligned_bam/{sample}_aln.bam", prefix=PREFIX, sample=SAMPLE),
        sorted_bam_reads= expand("results/{prefix}/{sample}/post_align/sorted_bam/{sample}_aln_sort.bam", prefix=PREFIX, sample=SAMPLE),
        bam_reads_duplicates_removed = expand("results/{prefix}/{sample}/post_align/remove_duplicates/{sample}_aln_marked.bam", prefix=PREFIX, sample=SAMPLE),
        dups_rmvd_sorted_bam_reads = expand("results/{prefix}/{sample}/post_align/sort_bam_dups_rmvd/{sample}_aln_duplicates_removed_sort.bam", prefix=PREFIX, sample=SAMPLE),
        alignment_stats = expand("results/{prefix}/{sample}/coverage_depth/{sample}_alignment_stats.tsv", prefix=PREFIX, sample=SAMPLE),
        gatk_DoC = expand("results/{prefix}/{sample}/coverage_depth/{sample}_depth_of_coverage.sample_summary", prefix=PREFIX, sample=SAMPLE),
        unmapped_bam = expand("results/{prefix}/{sample}/bedtools/bedtools_unmapped/{sample}_unmapped.bed", prefix=PREFIX, sample=SAMPLE),
        bioawk_ref_size_file = expand("results/{prefix}/bioawk/{ref_name}.size", prefix=PREFIX, ref_name=get_ref_name(config["reference_genome"]))
        #unmapped_bam_positions = expand("results/{prefix}/{sample}/bedtools/bedtools_unmapped/{sample}_unmapped.bed_positions", prefix=PREFIX, sample=SAMPLE)
        
# trims the raw fastq files to give trimmed fastq files
rule clean:
    input:
        r1 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
    output:
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
        threads = config["ncores"] # this uses bwa ncores
    log:
        trim_log = "logs/{prefix}/{sample}/trimmomatic/{sample}.log"
    conda:
        "envs/trimmomatic.yaml"
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
        aligned_sam_out = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam"
    params:
        #outdir = "results/{prefix}/{sample}/align_reads",
        num_cores = config["ncores"],
        ref_genome = config["reference_genome"]
        #prefix = "{sample}"
    log:
        bwa_log= "logs/{prefix}/{sample}/align_reads/{sample}.log"
    conda:
        "envs/bwa.yaml"
    #wrapper:
        # how to call bwa prepare readgroup python file here
    #run:
        #read_group = prepare_readgroup({input.r1})
        #shell("bwa mem -M -R {read_group} -t {params.num_cores} {params.ref_genome} {input.r1} {input.r2} > {output.aligned_sam_out}") 
    shell:
        """
        split_field=$(python3 -c "from python_scripts.prepare_readgroup import prepare_readgroup; print(prepare_readgroup('{input.r1}'))") &&
        bwa mem -M -R "$split_field" -t {params.num_cores} {params.ref_genome} {input.r1} {input.r2} > {output.aligned_sam_out}
        """
        
# samclip
rule post_align_samclip:
    input:
        aligned_sam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/align_reads/{wildcards.sample}_aln.sam")
    output:
        clipped_sam_out = f"results/{{prefix}}/{{sample}}/post_align/samclip/{{sample}}_clipped.sam",
    conda:
        "envs/samclip.yaml"
    params:
        #outdir = "results/{prefix}/{sample}/post_align/sort_bam",
        #prefix = "{sample}",
        ref_genome= config["reference_genome"]
    shell:
        "samclip --ref {params.ref_genome} --max 10 < {input.aligned_sam_out} > {output.clipped_sam_out}" 
        
# sam to bam, sort bam file
rule post_align_sorted_bam:
    input:
        clipped_sam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/samclip/{wildcards.sample}_clipped.sam")
    output:
        bam_out = f"results/{{prefix}}/{{sample}}/post_align/aligned_bam/{{sample}}_aln.bam",
        sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sorted_bam/{{sample}}_aln_sort.bam"
    conda:
        "envs/samtools.yaml"
    params:
        outdir_temp = "results/{prefix}/{sample}/post_align/sorted_bam/{sample}_aln_sort_temp",
        prefix = "{sample}",
        ref_genome= config["reference_genome"]
    shell:
        """
        samtools view -Sb {input.clipped_sam_out} > {output.bam_out} &&
        samtools sort {output.bam_out} -m 500M -@ 0 -o {output.sorted_bam_out} -T {params.outdir_temp}
        """
        
# remove duplicates
rule post_align_remove_dups:
    input:
        sorted_bam_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sorted_bam/{wildcards.sample}_aln_sort.bam")
    output:
        bam_duplicates_removed_out = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam"
    params:
        outdir = "results/{prefix}/{sample}/post_align/remove_duplicates",
        prefix = "{sample}"
    conda:
        "envs/picard.yaml"
    shell:
        "picard MarkDuplicates REMOVE_DUPLICATES=true INPUT={input.sorted_bam_out} OUTPUT={output.bam_duplicates_removed_out} METRICS_FILE={params.outdir}/{params.prefix}_markduplicates_metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"
        
# sort and index bam file with duplicates removed      
rule post_align_sort_index_rmvd_dups:
    input:
        bam_duplicates_removed_out = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/remove_duplicates/{wildcards.sample}_aln_marked.bam")
    output:
        dups_rmvd_sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sort_bam_dups_rmvd/{{sample}}_aln_duplicates_removed_sort.bam"
    params:
        outdir = "results/{prefix}/{sample}/post_align/sort_bam_dups_rmvd/",
        prefix = "{sample}"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools sort {input.bam_duplicates_removed_out} -m 500M -@ 0 -o {output.dups_rmvd_sorted_bam_out} -T {params.outdir}/{params.prefix}_aln_sort_temp &&
        samtools index {output.dups_rmvd_sorted_bam_out}
        """      
        
# determine statistics of file
rule stats:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sort_bam_dups_rmvd/{wildcards.sample}_aln_duplicates_removed_sort.bam")
    output:
        alignment_stats = f"results/{{prefix}}/{{sample}}/coverage_depth/{{sample}}_alignment_stats.tsv" 
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools flagstat {input.index_sorted_dups_rmvd_bam} > {output.alignment_stats}" 

# determine coverage of bam file       
rule coverage_depth:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sort_bam_dups_rmvd/{wildcards.sample}_aln_duplicates_removed_sort.bam")
    output:
        gatk_depthCoverage_summary = f"results/{{prefix}}/{{sample}}/coverage_depth/{{sample}}_depth_of_coverage.sample_summary"
    params:
        outdir = "results/{prefix}/{sample}/coverage_depth",
        ref_genome = config["reference_genome"], 
        prefix = "{sample}",
        intervals = config["bed_file"]
    conda:
        "envs/gatk.yaml"
    #wrapper:
        #"file:python_scripts/gatk"
    shell:
        "gatk DepthOfCoverage -R {params.ref_genome} -O {params.outdir}/{params.prefix}_depth_of_coverage -I {input.index_sorted_dups_rmvd_bam} --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {params.intervals}"

# created ref genome files 
rule bioawk:
    input:
        ref_name=get_ref_name(config["reference_genome"]),
        #ref_name = lambda wildcards: get_ref_name(config["reference_genome"]),
        ref_genome=config["reference_genome"]
    output:
        reference_size_file=f"results/{{prefix}}/ref_genome_files/{{ref_name}}.size"
    conda:
        "envs/bioawk.yaml"
    shell:
        "bioawk -c fastx '{{ print $name, length($seq) }}' < {input.ref_genome} > {output.reference_size_file}"

#create bed file from biowk
#rule created_bed_file:
    #input:
        # reference_size_file = output from rule bioawk
    #output:
        #reference_window_file = f"results/{{prefix}}/ref_genome_files/{{ref_name}}.bed"
    #conda:
        #"envs/bedtools.yaml"
    #shell:
        #"bedtools makewindows -g {reference_size_file} -w 1000 > {reference_dir}/{first_part}.bed"

# bedtools
rule bedtools:
    input:
        index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sort_bam_dups_rmvd/{wildcards.sample}_aln_duplicates_removed_sort.bam")
    output:
        unmapped_bed = f"results/{{prefix}}/{{sample}}/bedtools/bedtools_unmapped/{{sample}}_unmapped.bed"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools genomecov -ibam {input.index_sorted_dups_rmvd_bam} -bga | awk '$4==0' > {output.unmapped_bed}"

# this rule isnt working
# returns unmapped positions file    
#rule parse_bed_file:
    #input:
        #unmapped_bed = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/bedtools/bedtools_unmapped/{wildcards.sample}_unmapped.bed")
    #output:
        #unmapped_bam_positions = f"results/{{prefix}}/{{sample}}/bedtools/bedtools_unmapped/{{sample}}_unmapped.bed_positions"
    #run:
        #parse_bed_file({input.unmapped_bed})
        
# this rule is not working   
#rule bedgraph_cov:
    #input:
        #index_sorted_dups_rmvd_bam = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/post_align/sort_bam_dups_rmvd/{wildcards.sample}_aln_duplicates_removed_sort.bam")
    #output:
        #bedgraph_cov = f"results/{{prefix}}/{{sample}}/bedgraph_coverage/{{sample}}.bedcov"
    #params:
        #ref_genome = config["reference_genome"]
    #conda:
        #"envs/bedtools.yaml"
    #wrapper:
     
    
    
    
    
    
    
    


