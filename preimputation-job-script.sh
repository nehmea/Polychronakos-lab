#!/bin/bash
#SBATCH -J preimputation-pipeline
#SBATCH --mem=8G
#SBATCH -t 3:00:00
#SBATCH --account=def-cpolychr

#Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -p|--prefix) PREFIX="$2"; shift;;       # Output file name
          -str|--strand_file) STR_FILE="$2"; shift;;  # Strand file name
    esac
    shift
done

#give access to script file
chmod +x ./preimputation-script.sh

#run script
./preimputation-script.sh -p $PREFIX -str $STR_FILE