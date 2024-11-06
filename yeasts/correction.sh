#!/bin/bash

ref_name="${1}" #es : msh2_BYBY
path_exec="${2}"
ploidy="2"

base_dir=$(cd "mulo-gpw/pileup/common/${ref_name}" && pwd)

cd "${base_dir}"

all_dir=$(ls -d */ | cut -d "/" -f 1)

if [[ -n "$all_dir" ]]; then
 for dir in ${all_dir}; do
  files_name=$(ls "${dir}/common/"*".dat")
  count_name=$(echo "${files_name}" | wc -w)
  name_seq=($(echo "${files_name}" | cut -d "-" -f 1 | cut -d "/" -f 3))

  [[ -d "${dir}/correction" ]] || mkdir -p "${dir}/correction"

  if [[ "$count_name" -ge 2 ]]; then
   #printf "%s\n" "${name_seq[0]}"
   if [[ "${name_seq[0]}" == "${dir}" ]]; then
    ancestor="${name_seq[1]}" 
    endpoint="${name_seq[0]}"
   else
    ancestor="${name_seq[0]}" 
    endpoint="${name_seq[1]}"
   fi
   #printf "Ancestor: %s\n" $ancestor
   #printf "Endpoint: %s\n" $endpoint
   modal_ancestor=$(cat "${dir}/modal_coverage/${ancestor}-modalcov.dat" | head -n 1)
   std_ancestor=$(cat "${dir}/mean_coverage/${ancestor}-mc.dat" | tail -n 1)

   modal_endpoint=$(cat "${dir}/modal_coverage/${endpoint}-modalcov.dat" | head -n 1)
   std_endpoint=$(cat "${dir}/mean_coverage/${endpoint}-mc.dat" | tail -n 1)
   (
   "${path_exec}/correction_yeast" "${dir}/common/${ancestor}-${ref_name}_common_control_mpileup.dat" "${dir}/common/${endpoint}-${ref_name}_common_evolved_mpileup.dat" "${dir}/correction/" "${modal_ancestor}" "${std_ancestor}" "${modal_endpoint}" "${std_endpoint}" "${ploidy}"
   )
  fi

 done 
fi
