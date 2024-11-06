#!/bin/bash

out_dir="${1}"

cd "${out_dir}"

file=$(ls "list_files/"*"txt")

while IFS= read -r line
do
 #echo "${line}"
 (
 fasterq-dump --progress -e 4 --split-files "${line}"
 )

done < "$file"
