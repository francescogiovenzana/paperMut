# correction
Correction computes the number of bases required for the Bonferroni correction in the p-value calculation.

## Input
In the **correction/exec/** folder, there are two files: one to correct CRC (colorectal cancer) and the other to correct yeasts.  
They read the same types of input values.

1. the path of the ancestor with common lines,
2. the path of the endpoint with common lines,
3. the output path (without specyfing the name of the output file),
4. modal coverage of the ancestor,
5. standard dev of the ancestor,
6. modal coverage of the endpoint,
7. standard dev of the endpoint,
8. ploidy.

## Output

Output is a txt file with the number of bases required for the Bonferroni correction in the p-value  
calculation when the list of possible mutations is computed.  


This correction is useful when the probability of observing a rare event increases, such as reading a particular nucleotide.  
This occurs in our data due to the large number of reads in each dataset.

This correction should be applied for each ploidy value.

