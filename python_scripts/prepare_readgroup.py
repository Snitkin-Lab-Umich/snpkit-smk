# prepare_readgroup.py

import os
import re
import gzip

def prepare_readgroup(forward_read):
    samplename = os.path.basename(forward_read)
    if forward_read.endswith("fastq.gz"):
        with gzip.open(forward_read, 'rb') as output:
            firstLine = output.readline().decode("utf-8")
            if ":" in firstLine:
                split_field = re.split(r":", firstLine)
                id_name = split_field[1].rstrip()
            elif "/" in firstLine:
                split_field = re.split(r"/", firstLine)
                id_name = split_field[1].rstrip()
            else:
                id_name = "1"
            split_field = f"@RG\\tID:{id_name}\\tSM:{samplename}\\tLB:1\\tPL:Illumina"
            return split_field
    else:
        raise ValueError('Sequence read file extension not recognized.')
