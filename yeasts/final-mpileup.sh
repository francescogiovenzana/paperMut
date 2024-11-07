#!/bin/bash

ref_type="${1}"  #es : control or evolved
ref_name="${2}"  #es : msh2_BYBY 
path_exec="${3}" #exec path for parsing progrma

base_dir=$(cd "data/pileup/${ref_type}" && pwd)

out_dir="data/pileup/${ref_type}/${ref_name}"

[[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

cd "${base_dir}"

all_file=$(ls *".dat")

if [[ -n "$all_file" ]]; then
 for file in ${all_file}; do
  #printf "%s\n" $file
  name=$(echo "${file}" | cut -d "-" -f 1)
  name_fin="${name}-${ref_name}_${ref_type}_mpileup.dat"
  (
  "${path_exec}/parsing" "${file}" "${base_dir}/${ref_name}/${name_fin}"
  )
 done 
fi

#rm -f *".dat"

