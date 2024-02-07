shell("bioawk -c fastx '{ print $name, length($seq) }' < {snakemake.params.ref_genome} > {snakemake.output.reference_SIZE_file}")
