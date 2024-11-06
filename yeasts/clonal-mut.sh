#!/bin/bash

ref_name="${1}" #es : msh2_BYBY
path_exec="${2}"
n_gen="1000"
ploidy=2

base_dir=$(cd "mulo-gpw/pileup/common/${ref_name}" && pwd)

cd "${base_dir}"

all_dir=$(ls -d */ | cut -d "/" -f 1)

if [[ -n "$all_dir" ]]; then
 for dir in ${all_dir}; do
  files_name=$(ls "${dir}/common/"*".dat")
  count_name=$(echo "${files_name}" | wc -w)
  name_seq=($(echo "${files_name}" | cut -d "-" -f 1 | cut -d "/" -f 3))

  file_reeds="file_reeds_p${ploidy}.dat"
  file_no_mut="file_no_mutations_p${ploidy}.dat"

  #printf "%s\n" "${file_reeds}"
  #printf "%s\n" "${file_no_mut}"

  [[ -d "${dir}/reads/mr_clonal" ]] || mkdir -p "${dir}/reads/mr_clonal"

  if [[ "$count_name" -ge 2 ]]; then


   if [[ ! -e "${dir}/generations.dat" ]]; then
    touch "${dir}/generations.dat"
    echo "$n_gen" > "${dir}/generations.dat"
   fi

   n_gen_final=$(cat "${dir}/generations.dat" | head -n 1)
   #printf "%s\n" "${n_gen_final}"

   if [[ ${ploidy} -gt 1 ]]; then
    printf "%s\n" "${dir}"
    (
    "${path_exec}/clonal" "${dir}/reads/${file_reeds}" "${dir}/reads/${file_no_mut}" "${dir}/reads/mr_clonal/" "${n_gen_final}"
    )
   else

    if [[ "${name_seq[0]}" == "${dir}" ]]; then
     ancestor="${name_seq[1]}" 
     endpoint="${name_seq[0]}"
    else
     ancestor="${name_seq[0]}" 
     endpoint="${name_seq[1]}"
    fi

    modal_endpoint=$(cat "${dir}/modal_coverage/${endpoint}-modalcov.dat" | head -n 1)
    printf "%s\n" "hap:${dir}"
    printf "%s\n" "${modal_endpoint}"
    (
    "${path_exec}/clonal_hap" "${dir}/reads/${file_reeds}" "${dir}/reads/${file_no_mut}" "${dir}/reads/mr_clonal/" "${modal_endpoint}" "${n_gen_final}"
    )
   fi

  fi

 done 

fi
