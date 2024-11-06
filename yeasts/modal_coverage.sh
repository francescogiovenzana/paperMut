#!/bin/bash

ref_name="${1}" #es : msh2_BYBY
path_exec="${2}"
max_thr="1000"

base_dir=$(cd "mulo-gpw/pileup/common/${ref_name}" && pwd)

cd "${base_dir}"

all_dir=$(ls -d */ | cut -d "/" -f 1)
#
if [[ -n "$all_dir" ]]; then
 for dir in ${all_dir}; do
  files_name=$(ls "${dir}/common/"*".dat")
  printf "%s\n" $dir

  [[ -d "${dir}/modal_coverage" ]] || mkdir -p "${dir}/modal_coverage"

  if [[ -n "$files_name" ]]; then
   for f in ${files_name}; do
    #printf "%s\n" $f
    file_out="$(echo "${f}" | cut -d "-" -f 1 | cut -d "/" -f 3)-modalcov.dat"
    #printf "%s\n" $file_out
    (
    "${path_exec}/modal_coverage" "${f}" "${dir}/modal_coverage/${file_out}"
    )
   done
  fi

 done 
fi
