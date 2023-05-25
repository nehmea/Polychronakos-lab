#!/bin/bash
#Read the argument values
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -i|--input) input_name="$2"; shift;;
      -o|--output) output_name="$2"; shift;;
	-str|--strand_file) str_file="$2"; shift;;
    esac
    shift
done

#run this command 
plink --file $input_name --make-bed --out $input_name

#CHIP: GSAMD-24v3-0-EA_20034606_A1-b37.strand
#ORIENTATION: TOP strand to Forward strand
chr_file=GSAMD-24v3-0-EA_20034606_A1-b37.strand.chr
pos_file=GSAMD-24v3-0-EA_20034606_A1-b37.strand.pos
flip_file=GSAMD-24v3-0-EA_20034606_A1-b37.strand.flip
cat $str_file | cut -f 1,2 > GSAMD-24v3-0-EA_20034606_A1-b37.strand.chr
cat $str_file | cut -f 1,3 > GSAMD-24v3-0-EA_20034606_A1-b37.strand.pos
cat $str_file | awk '{if ($5=="-") print $0}' | cut -f 1 > GSAMD-24v3-0-EA_20034606_A1-b37.strand.flip


#Print the argument values
plink --allow-no-sex --bfile $input_name --update-map $chr_file --update-chr --make-bed --out TEMP_FILE_1
plink --allow-no-sex --bfile TEMP_FILE_1 --update-map $pos_file --make-bed --out TEMP_FILE_2
plink --allow-no-sex --bfile TEMP_FILE_2 --flip $flip_file --make-bed --out TEMP_FILE_3
plink --allow-no-sex --bfile TEMP_FILE_3 --extract $pos_file --make-bed --out $output_name
del -f TEMP_FILE_*


##create plink text files
#return to local
plink --bfile $output_name --recode --out $output_name


#recode plink text files into VCF
plink --file $output_name --recode vcf --out $output_name