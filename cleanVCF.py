# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:06:57 2022

This script removes specific chromosomes from a VCF file and saves the cleaned data in a new file.
Authors: alixv, nehmea
"""

import sys

# Get the filename from command-line argument
filename = sys.argv[1]  # prints python_script.py

# Open the output file for writing
output = open(f"{filename}_clean.vcf", "w")

# List of chromosomes to exclude
chros = ["0", "23", "24", "25", "26", "MT", "X", "XY", "Y"]

# Open the input VCF file for reading
with open(f"{filename}.vcf", "r") as f:
    print("START")

    # Iterate over each line in the input file
    for line in f:
        if line[0] != "#":
            # Split the line into data elements
            data = line.strip().split()

            # Check if the chromosome is in the list of excluded chromosomes
            if data[0] not in chros:
                # Write the line to the output file
                output.write(line)
        else:
            # Write the header lines to the output file
            output.write(line)

# Close the output file
output.close()

print("DONE")
