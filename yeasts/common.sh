#!/bin/bash

ref_name="${1}" #es : msh2_BYBY
path_exec="${2}"

base_dir=$(cd "data/pileup" && pwd)

out_dir="common/${ref_name}"

cd "${base_dir}"

all_file=$(ls "evolved/${ref_name}/"*".dat")
c_file=$(ls "control/${ref_name}/"*".dat")
c_name=$(echo "${c_file}" | cut -d "-" -f 1 | cut -d "/" -f 3)
name_a="${c_name}-${ref_name}_common_control_mpileup.dat"

#printf "%s\n" $c_name

if [ -n "$all_file" ] && [ -n "$c_file" ]; then
 for file in ${all_file}; do
  name=$(echo "${file}" | cut -d "-" -f 1 | cut -d "/" -f 3)
  #printf "%s\n" $c_file
  #printf "%s\n" $file
  [[ -d "${out_dir}/${name}/common" ]] || mkdir -p "${out_dir}/${name}/common"
  name_e="${name}-${ref_name}_common_evolved_mpileup.dat"
  #printf "%s\n" $name_a
  #printf "%s\n" $name_e
  (
  "${path_exec}/common" "${c_file}" "${file}" "${out_dir}/${name}/common/${name_a}" "${out_dir}/${name}/common/${name_e}"
  )
 done 
fi

