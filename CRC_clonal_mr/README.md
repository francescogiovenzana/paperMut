# CRC clonal mutation rate
This program computes the mutation rate by counting only clonal mutations,  
those mutations that, in a diploid organism, appear at a frequency around 0.5.

## Input

1. the path of ``file_reeds_p2.dat`` (possible mutations for CN = 2),
2. the path of ``file_no_mutations_p2.dat`` (number of non mutated bases for CN = 2),
3. the path of ``file_reeds_p3.dat`` (possible mutations for CN = 3),
4. the path of ``file_no_mutations_p3.dat`` (number of non mutated bases for CN = 3),
5. the path for the output files,
6. the clonal number (1307, 1502 or 0282 are the numbers that identify the different clones of the MAL lines),
7. the number that identify different lines of the same clone (9 for 1307, 3 for 1502, 5 for 0282).

Number (6) and (7) are mandatory to select the correct number of generations for each MAL experiments.

The mutation rate is then easily computed as

```math
\frac{\text{mutated bases}}{\text{total number of bases} * \text{generations}}
```

## Output

The output are 2 txt files:

1. ``file_mr_clonal.dat``
2. ``file_n_mutated_bases.dat``

File (1) contains the clonal mutation rates computed at different sigma levels,  
where sigma represents the standard deviation of a binomial distribution with

```math
\text{probability} = \frac{1}{CN}  
```

and

```math
\text{N} = \text{total coverage of a single base}.
```

The file includes mutation rate values for the following sigma levels (in this order):

```bash
0.1 sigma
0.3 sigma 
1 sigma
2 sigma
3 sigma
4 sigma
5 sigma
6 sigma
```
File (2) contains the number of selected mutations at each sigma level, in the same order.
