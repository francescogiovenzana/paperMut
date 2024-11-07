#!/bin/bash

ref_name="S288C_BY"
#n_ploidy="2"

## system's settings ----------------------------------------------------------

### check logs folder
if [[ ! -d "logs" ]]; then mkdir "logs"; fi

### short-reads subshell
(
### reference indexing
bash index-ref.sh "${ref_name}" > "logs/index-ref.out" 2> "logs/index-ref.err"

### mapping
bash map-sr.sh "${ref_name}" > "logs/map-sr.out" 2> "logs/map-sr.err"
)
