# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 16:06:57 2022

@author: alixv
@author2: nehmea
"""
import sys

filename = sys.argv[1]  # prints python_script.py

output = open(f"{filename}_clean.vcf", "w")
chros = ["0", "23", "24", "25", "26", "MT", "X", "XY", "Y"]

with open(f"{filename}.vcf", "r") as f:
    print("START")
    for line in f:
        if line[0] != "#":
            data = line.strip().split()
            if data[0] not in chros:
                output.write(line)
        else:
            output.write(line)

output.close()

print("DONE")
