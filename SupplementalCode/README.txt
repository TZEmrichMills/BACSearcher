Supplemental Code

Developed by Tom Emrich-Mills and John W. Davey

Associated manuscript: A recombineering pipeline to clone large and complex genes in Chlamydomonas, by Emrich-Mills et al. (2020)

------

FILES


BACSearcher.txt

The BACSearcher python script, including brief usage instructions.
Requires gffutils, biopython, primer3-py and IntervalTree modules.
For more detailed information regarding function, usage and precautions, see Supplemental Method 2.

-

UnsplicedGenes.txt

Python script for extracting the sequences, lengths and start/end positions of Chlamydomonas nuclear genes from the full genome fasta Creinhardtii_281_v5.fa.gz (see below).
Contains brief usage instructions. 

-

Creinhardtii_281_v5.0.fa.gz - NOT PROVIDED

Precursor file required by BACSearcher.txt and UnsplicedGenes.txt.
Contains a fasta file of the Chlamydomonas genome.
Can be downloaded from Phytozome>PhytozomeV12>Creinhardtii>assembly. 
Download may require login.

-

Creinhardtii_281_v5.5.gene.gff3.gz - NOT PROVIDED

Precursor file required by BACSearcher.txt and UnsplicedGenes.txt.
Contains annotation information for the 17741 Chlamydomonas nuclear genes downloaded from Phytozome 12.
Can be downloaded from Phytozome>PhytozomeV12>Creinhardtii>annotation. 
Download may require login.
-

BACs_fosmids.pairs.tsv

Precursor file required by BACSearcher to predict BAC coverage for each gene.
For usage instructions, see Supplemental Method 2A; this file corresponds to precursor file I.
The file contains start and end coordinates for all valid BACs and fosmids within the Chlamydomonas BAC library, as well as the lengths of each. 
BACs are identified by a PTQ number and fosmids are identified by a VTP number. 
Plasmids are defined as valid if their start and end sequences are mapped to the same chromosome, with one end mapped to the sense strand and the other to the anti-sense strand. 
These criteria exclude 3179 of the 11,676 BACs (27.2%) and 20,696 of the 56,276 fosmids (36.8%) that are annotated against the Chlamydomonas genome (v5.5). 
Due to the way the BACs have been annotated against the genome, some valid plasmids have more than two ends mapped to the same chromosome and so are represented more than once (see Supplemental Method 2B for additional information).

-

BAC_wells.txt

Precursor file required by BACSearcher in order to provide plate coordinates for any predicted BACs within the Chlamydomonas BAC library.
For usage instructions, see Supplemental Method 2A; this file corresponds to precursor file V. 

-

TRF-output-sorter.txt

Python script intended for processing the output information produced by Tandem Repeats Finder (Benson, 1999) into CSV format. 
Usage instructions are included at the top of the script and in Supplemental Method 4A.
See Supplemental Method 3A for details on Tandem Repeats Finder usage.

-

WindowMasker-output-sorter.txt

Python script intended for counting the number of repeat sequences masked by WindowMasker (Morgulis et al., 2006).
The script takes the interval TXT file produced by WindowMasker and outputs a CSV file containing the number or repeats per input sequence.
Usage instructions are included at the top of the script and in Supplemental Method 4B.
See Supplemental Method 3C for details on WindowMasker and its usage.

-

Check_primers-output-sorter.txt

Python script intended to enable quick analysis of the output from the Primer3 check_primers module (Rozen and Skaletsky, 2000). 
Usage instructions are included at the top of the script and in Supplemental Method 4C. 
See Supplemental Method 3D for details of the Primer3 check_primers module.

------

EXAMPLE DATA


Example_list.txt

Example list of four genes. Note the format of gene IDs required for BACSearcher and UnsplicedGenes.

-

Example_BACSearcher

BACSearcher output files generated using Example_list.txt. Includes .tsv files for suitable BACs (.bacs.tsv) and fosmids (.fosmids.tsv) for each gene.
Input: 
./BACSearcher.txt -p BACs_fosmids.pairs.tsv -f Creinhardtii_281_v5.0.fa.gz -g Creinhardtii_281_v5.5.gene.gff3.gz -l Example_list.txt -w Bac_wells.txt -o Example_BACSearcher

-

Example_UnsplicedGenes

UnsplicedGenes outputs generated using Example_list.txt. Includes .txt fasta files for 3'UTRs, 5'UTRs, ATG-Stop regions and full genes (FullSeq.txt). Also includes a tsv file showing positions and lengths (details.tsv).
Input: 
./UnsplicedGenes.txt -f Creinhardtii_281_v5.0.fa.gz -g Creinhardtii_281_v5.5.gene.gff3.gz -l Example_list.txt -o Example_UnsplicedGenes





