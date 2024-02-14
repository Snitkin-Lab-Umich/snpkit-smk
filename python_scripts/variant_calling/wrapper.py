__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell
#from os import path
#import re

ref_genome = snakemake.params.get("ref_genome", "")
mpileup_params = snakemake.params.get("mpileup_params", "")

#log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# what is happening here  
shell("bcftools mpileup -f {ref_genome} {snakemake.input.index_sorted_dups_rmvd_bam} | bcftools call -Ov -v -c -o {snakemake.output.final_raw_vcf}")
