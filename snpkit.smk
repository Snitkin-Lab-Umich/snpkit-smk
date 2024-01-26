# Author Ali Pirani and Dhatri Badri

import pandas as pd
import os

# How to read in R1 and R2 samples?
SAMPLES, READS = glob_wildcards('{sample}_{read_pair}_001.fastq.gz')
READS = list(set(READS)) # set is a list that consists of unique values 

#samples_df = pd.read_csv(config["samples"])
#SAMPLE = list(samples_df['sample_id'])
#PREFIX = config["prefix"]
#SHORTREADS = list(samples_df['sample_id'])

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

# Not sure if this function works/is correct
# downsample reads
def downsample(R1_file, R2_file, R1_out, R2_out, genome_size):
    
    R1_file = R1_file.pop()
    R2_file = R2_file.pop()
    R1_out = R1_out.pop()
    R2_out = R2_out.pop()

    # QUESTION: do you still wanna use mash to estimate genome size??
    # Run Mash to estimate Genome size
    #mash_cmd = "mash sketch -o /tmp/sketch_out -k 32 -m 3 -r %s >& /tmp/sketch_stdout" % R1_file
    #if genome_size:
        #gsize = int(genome_size)   
    #else:
        #try:
            #call(mash_cmd)
        #except sp.CalledProcessError:
            #sys.exit(1)

        #with open("/tmp/sketch_stdout", 'r') as file_open:
            #for line in file_open:
                #if line.startswith('Estimated genome size:'):
                    #gsize = float(line.split(': ')[1].strip())
                #if line.startswith('Estimated coverage:'):
                    #est_cov = float(line.split(': ')[1].strip())
        #file_open.close()

    gsize = genome_size.pop()

    print("Using Genome Size: %s to calculate coverage" % gsize)

    # Extract basic fastq reads stats with seqtk
    seqtk_check = "seqtk fqchk -q3 %s > /tmp/%s_fastqchk.txt" % (R1_file, R1_file)

    print(seqtk_check)
    try:
        os.system(seqtk_check)
    except sp.CalledProcessError: 
        print('Error running seqtk for extracting fastq statistics.')
        sys.exit(1)

    with open("%s_fastqchk.txt" % R1_file, 'r') as file_open: 
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
            print("seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R1_file, factor, r1_sub)) # what is seqtk sample path??
            seqtk_downsample = "seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R1_file, factor, r1_sub)) # what is seqtk sample path??
            os.system(seqtk_downsample)
        except sp.CalledProcessError:
            print('Error running seqtk for downsampling raw fastq reads.')
            sys.exit(1)

        if R2_file:
            #r2_sub = "/tmp/%s" % os.path.basename(R2_file)
            r2_sub = R2_out
            try:
                print("seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R2_file, factor, r2_sub))  # what is seqtk sample path??
                os.system("seqtk sample %s %s | pigz --fast -c -p 2 > %s" % (R2_file, factor, r2_sub))  # what is seqtk sample path??
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
    return r1_sub, r2_sub

# parse bed file 
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

rule all:
    trimmed_reads=expand('trimmomatic/{sample}_{reads}_trim_paired.fastq.gz', sample=SAMPLES, reads=READS),


# trims the raw fastq files to give trimmed fastq files
rule clean:
    input:
        r1 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["input_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz"))
    output:
        r1 = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_paired.fastq.gz",
        r2 = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_paired.fastq.gz", 
        r1_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R1_trim_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/{{sample}}/trimmomatic/{{sample}}_R2_trim_unpaired.fastq.gz",
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
        colon=config["colon"],
        target_length=config["targetlength"],
        crop_length=config["crop_length"]
        #threads = config["ncores"],
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
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        outr1 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R1_trim_paired.fastq.gz",
        outr2 = f"results/{{prefix}}/{{sample}}/downsample/{{sample}}_R2_trim_paired.fastq.gz",
    params:
        gsize = config["genome_size"]
    logs:
        downsample_log = "logs/{prefix}/{sample}/downsample/{sample}.log"
    run:
        downsample_reads({input.r1}, {input.r2}, {output.outr1}, {output.outr2}, {params.gsize}) &>{logs.downsample_log} 


# aligns trimmed fastq files using bwa to give sam file
rule align_reads:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R1_trim_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/{wildcards.sample}/trimmomatic/" + f"{wildcards.sample}_R2_trim_paired.fastq.gz"),
    output:
        aligned_sam_out = f"results/{{prefix}}/{{sample}}/align_reads/{{sample}}_aln.sam"
    params:
        outdir = "results/{prefix}/{sample}/align_reads",
        num_cores = config["ncores"],
        ref_genome = config["reference_genome"]
        #prefix = "{sample}",
    log:
        bwa_log= "logs/{prefix}/{sample}/align_reads/{sample}.log"
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem -M -R {input.R1} -t {params.num_cores} {params.ref_genome} {input.R1} > {output.aligned_sam_out}" # this command is not right 

# samclip
rule post_align_samclip:
    input:
        aligned_sam_out = lambda wildcards: expand(f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln.sam")
    output:
        clipped_sam_out = f"results/{{prefix}}/{{sample}}/post_align/samclip/{{sample}}_clipped.sam"
    conda:
        "envs/samtools.yaml"
    params:
        outdir = "results/{prefix}/{sample}/post_align/samclip"
        ref_genome= config["reference_genome"]
    shell:
        "samclip --ref {params.ref_genome} --max 10 < {input.aligned_sam_out} > {output.clipped_sam_out}"

#sam to bam
rule post_align_aln_bam:
    input:
        clipped_sam_out = f"results/{{prefix}}/{{sample}}/post_align/samclip/{{sample}}_clipped.sam"
    output:
        bam_out = f"results/{{prefix}}/{{sample}}/post_align/aligned_bam/{{sample}}_aln.bam"
    params:
        outdir = "results/{prefix}/{sample}/post_align/aligned_bam"
    conda:
        "envs/samtools.yaml" 
    shell:
        "samtools view -Sb {input.clipped_sam_out} > {output.bam_out}"  

# sort bam file
rule post_align_sort_bam:
    input:
        bam_out = f"results/{{prefix}}/{{sample}}/post_align/aligned_bam/{{sample}}_aln.bam"
    output:
        sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sorted_bam/{{sample}}_aln_sort.bam"
    params:
        outdir = "results/{prefix}/{sample}/post_align/sort_bam",
        #temp_outdir = "results/{prefix}/{sample}/post_align/sort_bam/{sample}_aln_sort_temp"
        prefix = "{sample}"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort {input.bam_out} -m 500M -@ 0 -o {output.sorted_bam_out} -T {params.outdir}/{params.prefix}_aln_sort_temp" 

# remove duplicates from sorted bam file
rule post_align_rm_duplicates:
    input:
        sorted_bam_out = lambda wildcards: expand(f"results/{{prefix}}/{{sample}}/post_align/sorted_bam/{{sample}}_aln_sort.bam")
    output:
        bam_duplicates_removed_out = f"results/{{prefix}}/{{sample}}/post_align/remove_duplicates/{{sample}}_aln_marked.bam"
    params:
        outdir = "results/{prefix}/{sample}/post_align/remove_duplicates",
        prefix = "{sample}
    conda:
        "envs/picard.yaml" # is this the correct way to download picard?
    shell:
        "picard MarkDuplicates REMOVE_DUPLICATES=true INPUT={input.sorted_bam_out} OUTPUT={output.bam_duplicates_removed_out} METRICS_FILE={params.outdir}/{params.prefix}_markduplicates_metrics CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT" 

# sort bam files that had duplicates removed
rule post_align_sort_removed_duplicates_bam:
    input:
        bam_duplicates_removed_out = lambda wildcards: expand(f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln_marked.bam") 
    output:
        dups_rmvd_sorted_bam_out = f"results/{{prefix}}/{{sample}}/post_align/sort_bam_dups_rmvd/{{sample}}_aln_duplicates_removed_sort.bam"
    params:
        outdir = "results/{prefix}/{sample}/post_align/sort_bam_dups_rmvd/",
        prefix = "{sample}"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort {input.bam_duplicates_removed_out} -m 500M -@ 0 -o {output.dups_rmvd_sorted_bam_out} -T {params.outdir}/{params.prefix}_aln_sort_temp"

#index bam file
rule post_align_index_bam:
    input:
        dups_rmvd_sorted_bam_out = lambda wildcards: expand(f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln_duplicates_removed_sort.bam")
    output:
        index_bam_out = f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln_duplicates_removed_sort_index.bam" # not referencing this in shell command?
    params:
        outdir = "results/{prefix}/{sample}/post_align/index_bam"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input.dups_rmvd_sorted_bam_out}"

# determine coverage and statistics of bam file
rule coverage_depth_stats:
    input:
        index_bam_out = lambda wildcards: expand(f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln_duplicates_removed_sort_index.bam")
    output:
        gatk_depthCoverage_summary = f"results/{{prefix}}/{{sample}}/coverage_depth/%s/%s_depth_of_coverage.sample_summary" % (out_path, analysis_name)
        alignment_stats = f"results/{{prefix}}/{{sample}}/coverage_depth/{{sample}}_alignment_stats" 
    params:
        outdir = "results/{prefix}/{sample}/coverage_depth",
        ref_genome= config["reference_genome"], # am i referencing the right ref genome here?
        prefix = "{sample}"
    conda:
        "envs/samtools.yaml",
        # how to download gatk? by zip or conda?
    shell:
        "samtools flagstat {input.index_bam_out} > {output.alignment_stats}" 
    run:
        interval = pattern.sub(lambda m: rep[re.escape(m.group(0))], {params.ref_genome})
        shell("gatk DepthOfCoverage -R {params.ref_genome} -O {params.outdir}/{params.prefix}_depth_of_coverage -I {input.index_bam_out} --summary-coverage-threshold 1 --summary-coverage-threshold 5 --summary-coverage-threshold 9 --summary-coverage-threshold 10 --summary-coverage-threshold 15 --summary-coverage-threshold 20 --summary-coverage-threshold 25 --ignore-deletion-sites --intervals {interval}")

rule bioawk:
    input:
        ref_genome = config["reference_genome"]
    output:
        bio_ref_size = "results/{prefix}/bioawk/"
    conda:
    shell:
        "bioawk -c fastx '{ print $name, length($seq) }' <{params.ref_genome} > {params.ref_genome}.size"

# this function is a little bit of a mess 
rule bedgraph_cov:
    input:
        index_bam_out = lambda wildcards: expand(f"results/{{prefix}}/{{sample}}/post_align/{{sample}}_aln_duplicates_removed_sort_index.bam")
    output:
        ref_size_file = 
        bed_out = f"results/{{prefix}}/{{sample}}/bedgraph_coverage/%s/%s.bed" % (reference_dir, first_part)
        bedgraph_coverage_out = ".bedcov"
        umapped_bed_file = f"results/{{prefix}}/{{sample}}/bedgraph_coverage/%s/%s_unmapped.bed" % (out_path, analysis)
        unmapped_positions_file = # _positions file?
    params:
        ref_genome= config["reference_genome"]
    conda:
    run:
        #shell("bioawk -c fastx '{ print $name, length($seq) }' <{params.ref_genome} > {params.ref_genome}.size")
        reference_SIZE_file = bioawk_make_reference_size({params.ref_genome})
        reference_filename_base = os.path.basename(reference)
        reference_first_part_split = reference_filename_base.split('.')
        first_part = reference_first_part_split[0]
        reference_dir = os.path.dirname({params.ref_genome})
        shell("bedtools makewindows -g {reference_SIZE_file} -w 1000 > {reference_dir}/{first_part}.bed")
        reference_windows_file = "%s/%s.bed" % (reference_dir, first_part)
        shell("bedtools coverage -abam {index.index_bam_out} -b {reference_windows_file} > {output.bedgraph_coverage_out}")
        shell("bedtools genomecov -ibam {input.index_bam_out} -bga | awk '$4==0' > {output.unmapped_bed_file}")
        final_bed_unmapped_file = {output.unmapped_bed_file}
        parse_bed_file(final_bed_unmapped_file) # what is the output from this file?



