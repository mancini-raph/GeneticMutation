## A2.py

This program compares SARS-CoV-2 reference gene sequences with variant gene sequences using data from a FASTA-formatted file. 

The script reads sequences from `sequencesfasta.txt`, separates reference (`NC`) sequences from variant (`ORF`) sequences, and compares matching genes codon-by-codon. 

For each gene, the program:

* Calculates the length of the reference and variant sequences
* Detects nucleotide differences between the sequences
* Determines whether each mutation is:

  * Silent (S): amino acid does not change
  * Missense (M): amino acid changes
  * Nonsense (N): mutation creates a stop codon
* Stores the first four detected mutations
* Counts the total number of sequence differences

Mutations are displayed in the format:
`A139G_M`

This means the nucleotide `A` at position `139` changed to `G`, producing a missense mutation.

The program uses a codon-to-amino-acid dictionary to translate codons and determine the biological effect of mutations. Finally, it outputs a formatted table containing the gene name, sequence lengths, detected mutations, and total mutation count.

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## a4.py

This program calculates evolutionary distances between DNA sequences using the Jukes–Cantor evolutionary model. 

The script reads DNA sequences from `sequence.txt` in FASTA format.  It then compares every sequence against all other sequences to determine how genetically similar or different they are.

The program:

* Parses sequence headers and DNA data
* Counts nucleotide differences between sequence pairs
* Normalizes the differences by sequence length
* Applies the Jukes–Cantor distance formula to estimate evolutionary distance
* Generates a symmetric distance matrix containing the distances between all species

The Jukes–Cantor model corrects for the possibility of multiple mutations occurring at the same nucleotide position over evolutionary time.

Smaller distance values indicate organisms that are more closely related genetically, while larger values indicate organisms that are more evolutionarily distant.

The final output is a formatted distance matrix showing the evolutionary relationship between all input sequences.
