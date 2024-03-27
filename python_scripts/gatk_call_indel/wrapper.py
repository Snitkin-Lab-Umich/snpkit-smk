__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell
#from os import path
#import re

ref_genome = snakemake.params.get("ref_genome", "")

#log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# 
shell("""
        gatk HaplotypeCaller -R {ref_genome} -I {snakemake.input.index_sorted_dups_rmvd_bam} -O {snakemake.output.final_raw_vcf} --native-pair-hmm-threads 8 &&
        gatk SelectVariants -R {ref_genome} -V {snakemake.output.final_raw_vcf} -select-type INDEL -O {snakemake.output.indel_file}
        """)
        
shell("""
       bgzip -c {snakemake.output.indel_file} > {snakemake.output.zipped_indel_vcf} 
       tabix -p vcf -f {snakemake.output.zipped_indel_vcf}
    """)