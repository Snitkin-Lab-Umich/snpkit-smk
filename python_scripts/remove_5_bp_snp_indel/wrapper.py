__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

import os
from snakemake.shell import shell
from os import path
#import re

raw_snp_vcf_file = os.path.basename(snakemake.input.snp_vcf[0])
raw_indel_vcf_file = os.path.basename(snakemake.input.indel_vcf[0])
output_file = os.path.basename(snakemake.output.remove_snps_5_bp_snp_indel_file)

print(raw_snp_vcf_file)
print(raw_indel_vcf_file)   
print(output_file)

def remove_5_bp_snp_indel(raw_snp_vcf_file, raw_indel_vcf_file):
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
    #print(type(raw_indel_vcf_file))
          
    return remove_snps_5_bp_snp_indel_file_name

#
#shell("""
       #bgzip -c {snakemake.output.remove_snps_5_bp_snp_indel_file} > {snakemake.output.zipped_indel_vcf} &&
       #tabix -p vcf -f {snakemake.output.zipped_indel_vcf}
    #""")