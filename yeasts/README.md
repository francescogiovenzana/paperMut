# yeasts
This folder contains various scripts useful for data analysis of yeast genome sequencing.

## Description
Except for the scripts that call some programs for indexing and mapping yeast genome sequences,  
the other scripts call the same programs found in the main folder and perform the same procedures  
explained in the individual folders, for all yeast ancestor and endpoint sequencing.

- ``index-ref.sh``: make indexes for ``bwa`` (script calls `bwa`).

Reference genome for Saccharomyces Cerevisiae used: ``S288C_reference_genome_R64-1-1_20110203``.

- ``map-sr.sh``: maps short-reads and fixes the bam file (script calls `bwa` and `samtools`).

- ``ref-map.sh``: calls both ``index-ref.sh`` and ``map-sr.sh`` in this order. It also creates  
``logs`` directory for outputs and errors of the used scripts.

- ``bam_to_pileup.sh``: transforms bam files in pileup files with `samtools mpileup`.

- ``final-mpileup.sh``: calls `parsing` to create all the pseudopileup txt files.

- ``common.sh``: calls `common` to create common lines files of the ancestor and endpoint sequencing.

- ``mean_coverage.sh``: calls `mean_coverage` to compute mean coverage and std dev for each yeast genome sequencing.

- ``modal_coverage.sh``: calls `modal_coverage` to compute modal coverage and for each yeast genome sequencing.

- ``correction.sh``: calls `correction_yeast` to compute the number of bases for the Bonferroni correction for every yeast genome sequencing.

- ``reads-mut.sh``: calls `list_mutations_yeast` to produce the list of possible mutations in each pair of ancestor-endpoint sequences.

- ``clonal-mut.sh``: calls `clonal` or `clonal_hap` to compute clonal mutation rate for each pair of ancestor-endpoint sequences.



### Notes

To maintain the folder structure as used, I have inserted a `.gitkeep` file in all the empty  
folders to upload them to GitHub. Before running any script, please remove these files from the folders.
