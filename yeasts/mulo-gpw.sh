#!/bin/bash

##concatenate multiple fastq files
#cat seq-a-1.fq seq-b-1.fq > allread-1.fq
#cat seq-a-2.fq seq-b-2.fq > allread-2.fq

## header ---------------------------------------------------------------------
##download fastq files with --progress and number of threads = 8
#fasterq-dump --progress -e 8 --split-files SRR12512594
#à

##transform bam in pileup format 
#samtools mpileup --min-MQ 5 --ignore-RG -f /path/to/N17.genome.fa /path/to/SRR9317924.bam > /path/to/output/SRR9317924.pileup.dat
##

### this script is the mulo-gpw runner

## user's settings ------------------------------------------------------------

ref_name="S288C_BY"
#n_ploidy="2"

### sample IDs where:
### - plus_samp is the (plus) selected sample in experimental design A or
###   the selected sample (plus or minus) in experimental design BH and BL
### - minu_samp is the (minus) selected sample in experimental design A or
###   the non-selected sample in experimental design BH and BL

### pm-1 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D2 D5 D8 D11"
### pm-2 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D5 D8 D11 D2"
### pm-3 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D8 D11 D2 D5"
### pm-4 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D11 D2 D5 D8"
### pe-1 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D3 D6 D9 D12"
### pe-2 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D6 D9 D12 D3"
### pe-3 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D9 D12 D3 D6"
### pe-4 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D12 D3 D6 D9"
### me-1 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D3 D6 D9 D12"
### me-2 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D6 D9 D12 D3"
### me-3 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D9 D12 D3 D6"
### me-4 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D12 D3 D6 D9"
### pz-1 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D13 D13 D13 D13"
### mz-1 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D13 D13 D13 D13"
### ez-1 (with exp_design="BH", assuming "equal" samples have a plus phenotype)
#plus_samp="D3 D6 D9 D12"
#minu_samp="D13 D13 D13 D13"

### the experimental design ("exp_design" variable) can be "A"
### (for plus/minus phenotyped samples)
### or "BH" (for plus-evolved/non-evolved samples)
### or "BL" (for minus-evolved/non-evolved samples)
#exp_design="BH"

## system's settings ----------------------------------------------------------

### check logs folder
if [[ ! -d "logs" ]]; then mkdir "logs"; fi

### short-reads subshell
(
### reference indexing
#bash index-ref.sh "${ref_name}" > "logs/index-ref.out" 2> "logs/index-ref.err"

### mapping
bash map-sr.sh "${ref_name}" > "logs/map-sr.out" 2> "logs/map-sr.err"

#bash gem.sh "${ref_name}" 150 > "logs/Mappability-GEM.out" 2> "logs/Time-Mappability-GEM.err"
)
