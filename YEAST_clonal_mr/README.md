# YEAST clonal mutation rate
This program computes the mutation rate by counting only clonal mutations,  
those mutations that, in a diploid organism, appear at a frequency around 0.5.

## Input

In **YEAST_clonal_mr/exec/** there are two files, ``clonal`` and ``clonal_hap`` where  
we splitted the case of diploid genomes from the case of haploid genomes.

Input parameters are basically the same:

1. the path of ``file_reeds_pX.dat`` (where X = copy number),
2. the path of ``file_no_mutations_pX.dat``,
3. the path for the output files,
4. modal coverage (only for ``clonal_hap``),
5. number of generations.


## Output 

Output files are exactly the same of those seen in **CRC clonal mutation rate**.


### Notes

We split the haploid case from the diploid case because we couldn’t calculate  
the error using a binomial distribution. We therefore used an alternative method to account for  
the selection error on the correct bases, taking the maximum between sequencing error and sensitivity error on the modal coverage.
