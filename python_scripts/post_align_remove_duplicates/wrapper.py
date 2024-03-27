__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell
#from os import path
#import re

outdir = snakemake.params.get("outdir", "")
outdir_dups_removed = snakemake.params.get("outdir_dups_removed", "")
prefix = snakemake.params.get("prefix", "") 
#ref_genome = snakemake.params.get("ref_genome", "")

#log = snakemake.log_fmt_shell(stdout=True, stderr=True)

#picard command - updated
shell("picard MarkDuplicates -REMOVE_DUPLICATES true -INPUT {snakemake.input.sorted_bam_out} -OUTPUT {snakemake.output.bam_duplicates_removed_out} -METRICS_FILE {outdir_dups_removed}/{prefix}_markduplicates_metrics -CREATE_INDEX true -VALIDATION_STRINGENCY LENIENT")

# samtools command
shell("samtools sort {snakemake.output.bam_duplicates_removed_out} -m 500M -@ 0 -o {snakemake.output.dups_rmvd_sorted_bam_out} -T {outdir}/{prefix}_aln_sort_temp"
      " && samtools index {snakemake.output.dups_rmvd_sorted_bam_out}")
