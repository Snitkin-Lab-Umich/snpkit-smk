__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
import re
from snakemake.shell import shell
#from snakemake_wrapper_utils.java import get_java_opts
from os import path


outdir = snakemake.params.get("outdir", "")
ref_genome = snakemake.params.get("ref_genome", "")
prefix = snakemake.params.get("prefix", "") 
#log = snakemake.log_fmt_shell(stdout=True, stderr=True)



#######################################

#shell("bioawk -c fastx '{ print $name, length($seq) }' < {snakemake.params.ref_genome} > {snakemake.output.reference_SIZE_file}")

reference_filename_base = os.path.basename({snakemake.input.ref_genome})
reference_first_part_split = reference_filename_base.split('.')
first_part = reference_first_part_split[0]
reference_dir = os.path.dirname({snakemake.input.ref_genome})

shell("bedtools makewindows -g {snakemake.output.reference_SIZE_file} -w 1000 > {reference_dir}/{first_part}.bed")

reference_windows_file = "%s/%s.bed" % (reference_dir, first_part)

shell("bedtools coverage -abam {snakemake.input.index_sorted_dups_rmvd_bam} -b {reference_windows_file} > %s/%s.bedcov" % (out_sorted_bam, reference_windows_file, out_path, analysis))









