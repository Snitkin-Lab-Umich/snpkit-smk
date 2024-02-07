__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts
from os import path

java_opts = get_java_opts(snakemake)
#extra = snakemake.params.get("extra", "")
outdir = snakemake.params.get("outdir", "")
ref_genome = snakemake.params.get("ref_genome", "")
prefix = snakemake.params.get("prefix", "") 
#log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Extract basename from the output file names
#out_basename = path.commonprefix(snakemake.output).rstrip(".")

rep = {"fasta": "bed", "fna": "bed"}
rep = dict((re.escape(k), v) for k, v in rep.items())
pattern = re.compile("|".join(rep.keys()))
interval = pattern.sub(lambda m: rep[re.escape(m.group(0))], {snakemake.params.ref_genome})

shell("gatk --java-options '{java_opts}' DepthOfCoverage" 
    "-R {snakemake.params.ref_genome}" 
    "-O {snakemake.params.outdir}/{snakemake.params.prefix}_depth_of_coverage" 
    "-I {snakemake.input.index_sorted_dups_rmvd_bam}"
    " --summary-coverage-threshold 1"
    " --summary-coverage-threshold 5" 
    "--summary-coverage-threshold 9" 
    "--summary-coverage-threshold 10" 
    "--summary-coverage-threshold 15" 
    "--summary-coverage-threshold 20" 
    "--summary-coverage-threshold 25" 
    "--ignore-deletion-sites"
    " --intervals {interval}"
)

