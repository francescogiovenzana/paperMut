# signatures
In the **signatures** directory, two files are implemented: one to extract the list of mutations with coordinates, and another to extract the list of possible signatures. 

## Input for the list of mutations with coordinates
In the **signatures/exec/** folder, there is one file to extract the list of possible mutations of CRC and yeast genomes with coordinates:  
- **list_coord_mut**

The input values are:

1. the path of the ancestor with common lines,
2. the path of the endpoint with common lines,
3. the output path (without specyfing the name of the output file),
4. modal coverage of the ancestor,
5. standard dev of the ancestor,
6. modal coverage of the endpoint,
7. standard dev of the endpoint,
8. copy number.

## Output list with coordinates
The output files are exactly the same as for the **list_mutations**, except for the coordinates of the possible mutations.

Here an example of output file:

```bash
..    ...       .   .   ... .   .
..    ...       .   .   ... .   .
14  60243132    A   C   103 1   2
14  60243134    T   G   106 1   2
14  60243141    A   T   104 1   2
14  60243145    T   A   100 1   2
14  60243160    T   G   103 1   2
14  60243162    T   G   107 1   2
14  60243164    G   A   112 1   2
14  60243167    T   G   108 1   2
14  60243171    T   A   112 2   2
14  60243178    T   G   110 1   2
14  60243199    A   C   118 1   2
14  60243292    T   G   97  1   2
14  60243374    G   A   106 1   2
14  60243382    A   C   112 1   2
14  60243421    A   C   107 1   2
14  60243422    G   T   102 1   2
14  60243435    A   C   107 1   2
14  60243443    A   C   117 1   2
14  60243461    C   G   100 1   2
14  60243505    A   C   110 1   2
14  60243546    T   C   100 1   2
14  60243689    G   A   85  1   2
14  60243881    A   T   90  1   2
14  60243982    A   C   88  1   2
14  60244027    T   G   89  1   2
14  60244051    C   A   93  1   2
14  60244053    G   T   94  1   2
14  60244095    C   A   109 1   2
..    ...       .   .   ... .   .
..    ...       .   .   ... .   .
```
where in the different columns you can find

- column 1) chromosome number
- column 2) base numberis the number of mutations 
- column 3) nucleotide on the ancestor
- column 4) nucleotide on the endpoint
- column 5) endpoint coverage
- column 6) number of mutated reads
- column 7) copy number

## Input for the signatures
In the **signatures/exec/** folder, there is one file to extract the list of the possible signatures:  
- **signatures_mut**

The input parameters are:

1. the path of the list of mutations with coordinates,
2. the path of the ancestor with common lines,
3. the path of the endpoint with common lines,
4. the output path (without specyfing the name of the output file),
5. copy number.

## Output signatures
An example of output file for the signatures:

```bash
.      ...      .   .   .   .   .   .   .   .   .
.      ...      .   .   .   .   .   .   .   .   .
4   113546389   C   C   G   C   T   G   C   C   G
4   113546398   G   C   T   G   A   T   G   C   T
4   113546423   C   C   A   C   A   A   C   C   A
4   113546435   T   T   T   T   A   T   T   T   T
4   113546501   T   G   C   T   T   C   T   G   C
4   113546557   T   A   A   T   G   A   T   A   A
4   113546653   T   T   T   T   A   T   T   T   T
4   113546723   A   A   A   A   T   A   A   A   A
4   113546733   A   C   A   A   G   A   A   C   A
4   113546736   G   T   G   G   G   TG  G   T   G
4   113546737   T   G   C   GT  T   C   T   G   C
4   113546743   C   A   T   C   T   T   C   A   T
4   113546746   T   T   A   T   A   A   T   T   A
4   113546778   T   C   G   T   G   G   T   C   G
4   113546793   C   A   C   C   T   C   C   A   C
4   113546812   C   A   C   C   T   C   C   A   C
4   113546826   G   A   G   G   T   TG  G   A   G
4   113546827   A   G   T   TA  T   T   A   G   T
4   113546838   A   T   G   A   G   G   A   T   G
4   113546877   A   G   C   A   T   C   A   G   C
4   113546895   T   G   A   T   T   A   T   G   A
4   113547676   T   G   G   T   A   G   T   G   G
4   113547811   C   A   G   C   G   G   C   A   G
4   113547932   T   T   C   T   C   C   T   T   C
4   113548168   A   T   AAT A   A   T   A   T   T
4   113548179   A   T   C   A   G   C   A   T   C
4   113548235   T   T   C   T   C   C   T   T   C
.      ...      .   .   .   .   .   .   .   .   .
.      ...      .   .   .   .   .   .   .   .   .
```
where in the different columns we can find

- column 1)  chromosome number
- column 2)  base number
- column 3)  previous ancestor nucleotide 
- column 4)  current ancestor nucleotide
- column 5)  next ancestor nucleotide
- column 6)  previous endpoint nucleotide
- column 7)  current endpoint nucleotide
- column 8)  next endpoint nucleotide
- column 9)  previous reference nucleotide
- column 10) current reference
- column 11) next reference nucleotide

Fundamentally, the analysis is performed on nucleotide triplets.
By "current," we mean the reference nucleotide where the mutation has been identified and to which the chromosome and base coordinates indicated here correspond. "Previous" and "next" refer to that nucleotide. The "reference" refers to the reads on the reference genome for those coordinates (hg38 for human and S288C for yeasts).


We emphasize again that, for the mutation on the current nucleotide, the one with the highest frequency was selected.

We have also reported when multiple occurrences appear in previous and next nucleotide for ancestor and endpoint. In the example above, it can be seen that for a given base, the nucleotide may not be unique, as the reads may have returned different nucleotides. These are highlighted in descending order of frequency: first, the one that was read most of the times, followed by the others that have at least one non-zero read.
e.g: in chr 4 base 113548168 the next on the ancestor could be A or T, but A has a higher frequency.
