#!/bin/bash

ref_name="${1}"

base_dir=$(cd "mulo-gpw/pileup/common/${ref_name}" && pwd)

cd "${base_dir}"

all_dir=$(ls -d */ | cut -d "/" -f 1)

for dir in ${all_dir}; do
 echo "$dir" >> "list_${ref_name}.txt"
done
