#!/bin/bash

##transform bam in pileup format 
##
seq_type="${1}" #es : control or evolved
mut_type="${2}" #es : msh2_BYBY 

base_dir=$(cd "mulo-gpw" && pwd)

cd "${base_dir}"

ref_genome="ref/S288C_BY-genome.fa"

all_fqs=$(ls "map-sr/"*".bam")

for f in ${all_fqs}; do
  #printf "%s\n" $f
  f_name=$(echo "${f}" | cut -d "/" -f 2 | cut -d "-" -f 1)
  #printf "%s\n" $f_name
  (
  samtools mpileup --min-MQ 5 --ignore-RG -f "${ref_genome}" "${f}" > "pileup/${seq_type}/${f_name}-${mut_type}_${seq_type}_pileup.dat"
  )
done

#rm -r "map-sr/"*".bam"
#rm -r "map-sr/"*".bam.bai"

