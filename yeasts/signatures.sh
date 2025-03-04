#!/bin/bash

ref_name="${1}"   #es : msh2_BYBY
path_exec="${2}"  #absolute path
ploidy="${3}"     #ploidy of the genome

printf "Ploidy = %s\n" $ploidy

base_dir=$(cd "data/pileup/common/${ref_name}" && pwd)

cd "${base_dir}"

all_dir=$(ls -d */ | cut -d "/" -f 1)

if [[ -n "$all_dir" ]]; then
 for dir in ${all_dir}; do
  files_name=$(ls "${dir}/common/"*".dat")
  count_name=$(echo "${files_name}" | wc -w)
  name_seq=($(echo "${files_name}" | cut -d "-" -f 1 | cut -d "/" -f 3))

  if [[ -d "${dir}/sign" ]]; then rm -rf "${dir}/sign"; fi
  mkdir -p "${dir}/sign"

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

   (
   "${path_exec}/signatures_mut" "${dir}/coord/file_reeds_p${ploidy}.dat" "${dir}/common/${ancestor}-${ref_name}_common_control_mpileup.dat" "${dir}/common/${endpoint}-${ref_name}_common_evolved_mpileup.dat" "${dir}/sign/" "${ploidy}"
   )
  fi

 done
fi
