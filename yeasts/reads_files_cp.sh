#!/bin/bash

ref_name="${1}"       #es : msh2_BYBY
dest_path="${2}"      #     destination directory path
directory_cp="${3}"   #     directory copied 

base_dir=$(cd "mulo-gpw/pileup/common/${ref_name}" && pwd)

cd "${base_dir}"

all_dir=$(ls -d */ | cut -d "/" -f 1)

if [[ -n "$all_dir" ]]; then
 for dir in ${all_dir}; do

  [[ -d "${dest_path}/${ref_name}/${dir}" ]] || mkdir -p "${dest_path}/${ref_name}/${dir}"

  cp -r "${dir}/${directory_cp}" "${dest_path}/${ref_name}/${dir}"

  printf "cp dir %s in dir %s\n" "${directory_cp}" "${dir}"

 done
fi
