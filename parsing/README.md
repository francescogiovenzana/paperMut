# parsing
Parsing converts a pileup-formatted text file into a more readable text format that also occupies less disk memory.


## Input
In the **parsing/exec/** folder, there are two files: one for parsing CRC (colorectal cancer) files and the other for parsing yeast files.

Input ``parsing``: 

1. the path of the input file pileup.gz of the CRC
2. the path of the output file pseudopileup.dat of the CRC

Input ``parsing_yeast``: 

1. the path of the input file pileup.dat of the yeast genome
2. the path of the output file pseudopileup.dat of the yeast genome

- Input file formats (example with a line in a CRC file):

``
chr1	12933843	12933844	G	31	...T..t.,tTTT.,.TT,,.t.,..T..t^].
``

where

chr1: chromosome number

12933843: 0-base

12933844: 1-base

G: reference

...T..t.,tTTT.,.TT,,.t.,..T..t^]. : forward and reverse reads.

Input file format for a yeast sequencing is the same, except for chromosome number: in the yeast  
sequencing chromosome are indicated with roman numerals. For e.g chr16 in yeast is XVI.

## Output

Output files are always txt files.

- Output file format (example with a line in a CRC file):

``
1	12933844	G	0	15	0	7	0	5	0	4
``

where we can see chr number, 1-base (you can select also 0-base, as you prefer), reference,  
reads forward (A,G,C,T) and reads reverse (a,g,c,t).


This procedure must be carried out for all ``ancestor`` and ``endpoint`` files of a MAL experiment, as well as any files divided by different ``CN (ploidy)``.
