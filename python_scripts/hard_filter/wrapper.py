__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell
#from os import path
#import re

#outdir = snakemake.params.get("outdir", "")
#outdir_dups_removed = snakemake.params.get("outdir_dups_removed", "")
#prefix = snakemake.params.get("prefix", "") 
#ref_genome = snakemake.params.get("ref_genome", "")
# params
ref_genome = snakemake.params.get("ref_genome", "")
dp_filter = snakemake.params.get("dp", "")
fq_filter = snakemake.params.get("fq", "")
mq_filter = snakemake.params.get("mq", "")
qual_filter = snakemake.params.get("qual", "")
af_filter = snakemake.params.get("af", "")

#log = snakemake.log_fmt_shell(stdout=True, stderr=True)

gatk_filter2_parameter_expression = "%s && %s && %s && %s && %s && %s" % (fq_filter, mq_filter, qual_filter, dp_filter, af_filter)
        
# snp 
shell("gatk VariantFiltration -R {ref_genome} -O {snakemake.output.filter2_snp_vcf} --variant {snakemake.input.final_raw_snp_vcf} --filter-expression \"{gatk_filter2_parameter_expression}\" --filter-name PASS_filter2")
shell("grep '#\|PASS_filter2' {snakemake.output.filter2_snp_vcf} > {snakemake.output.filter2_snp_final}")

# indels
shell("gatk VariantFiltration -R {ref_genome} -O {snakemake.output.filter2_indel_vcf} --variant {snakemake.input.final_raw_indel_vcf} --filter-expression \"{gatk_filter2_parameter_expression}\" --filter-name PASS_filter2")
shell("grep '#\|PASS_filter2' {snakemake.output.filter2_indel_vcf} > {snakemake.output.filter2_indel_final}")