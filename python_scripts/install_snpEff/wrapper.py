__author__ = "Dhatri Badri"
__copyright__ = "Copyright 2024, Dhatri Badri"
__email__ = "dhatrib@umich.edu"
__license__ = "MIT"

from snakemake.shell import shell
from os import path
import os
import re
import sys
import ruamel.yaml
from ruamel.yaml import YAML

# Directory path where .yaml files are located (relative path)
yaml_dir = os.path.join(".", ".snakemake", "conda")

# Initialize hash_name as None
hash_name = None

base_dir = snakemake.params.get("base_dir", "")
snpEff_db = snakemake.params.get("snpEff_db", "")

# Iterate through .yaml files
for filename in os.listdir(yaml_dir):
    if filename.endswith(".yaml"):
        yaml_path = os.path.join(yaml_dir, filename)
        #print(yaml_path)
        with open(yaml_path, 'r') as yaml_file:
            yaml_content = yaml_file.read()
            #print(yaml_content)
            # Search for "snpEff" in the file name
            if "snpEff" in yaml_content:
                # Extract the hash from the file name
                match = re.search(r'(\w+)\.yaml', filename)
                #print(match)
                if match:
                    hash_name = match.group(1)
                    #snpeff_path = os.path.join(yaml_dir, hash_name, "share", "snpeff-5.0-1", "snpEff.jar")
                    #print(f"Path of snpEff: {yaml_dir}/{hash_name}")
                    break  # Stop searching after finding the first occurrence

if 'hash_name' not in locals():
    print("snpEff hash name not found in .yaml files.")

snpEff_config_file = os.path.join(yaml_dir, hash_name, "share", "snpeff-5.0-1", "snpEff.config") 
snpeff_path = os.path.join(yaml_dir, hash_name, "share", "snpeff-5.0-1", "snpEff.jar")
data_directory = os.path.join(snakemake.params.get("base_dir", ""), "data")

shell(f"echo '{snakemake.params.snpEff_db}.genome : {snakemake.params.snpEff_db}' >> {snpEff_config_file}")
shell("java -jar {snpeff_path} build -genbank -v {snakemake.params.snpEff_db} -c {snpEff_config_file} -dataDir {data_directory}")


# Load the existing YAML content
yaml = ruamel.yaml.YAML(pure=True)
yaml.width = 100
# Set indentation parameters
#yaml.indent(mapping=0, sequence=0, offset=0)

with open("config/config.yaml", "r") as file:
    config = yaml.load(file)

# Update the specific values
config["snpEff_path"] = snpeff_path
config["data_dir"] = data_directory
config["snpEff_config_file"] = snpEff_config_file

# Write the updated content back to the file
with open("config/config.yaml", "w") as file:
    yaml.dump(config, file)