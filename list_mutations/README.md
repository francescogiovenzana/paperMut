# list mutations
List mutations returns the list of the possible mutations for a particular CRC or yeast genome.

## Input
In the **list_mutations/exec/** folder, there is one file to extract the list of possible mutations of CRC and yeast genomes:  
- **list_mutations**

The input values are:

1. the path of the ancestor with common lines,
2. the path of the endpoint with common lines,
3. the output path (without specyfing the name of the output file),
4. modal coverage of the ancestor,
5. standard dev of the ancestor,
6. modal coverage of the endpoint,
7. standard dev of the endpoint,
8. copy number.

## Output

The output are two txt files.  

The first txt file is called ``file_reeds_pX.dat``, with ``X`` equals the ploidy of the dataset.
This file contains the list of possible mutations and the data to compute the final frequency spectrum:

```bash
70      1       2
81      1       2
72      1       2
77      1       2
75      1       2
81      1       2
61      1       2
70      1       2
69      1       2
69      1       2
79      1       2
76      1       2
77      1       2
72      1       2
81      1       2
87      1       2
66      1       2
64      1       2
69      1       2
65      1       2
73      1       2
81      1       2
86      1       2
87      1       2
...    ...     ...
...    ...     ...
```
where

- first column is the total coverage
- second column is the number of mutations 
- third column is the copy number

The second txt file is called ``file_no_mutations_pX.dat``, with ``X`` equals the ploidy of the dataset.
This file contains the final number of non-mutated bases.

These files are present for each ploidy value in which the original dataset was divided (before performing any of the analyses described in the main folder).


