This repository contains two sets of genome analysis tools for the manuscript "A recombineering pipeline to clone large and complex genes in Chlamydomonas". The first is BACSearcher, a tool to facilitate easy recombineering of Chlamydomonas genes from the Chlamydomonas BAC library. The repository contains the BACSearcher python code plus two precursor files required run the script. The second set contains short python scripts that can be used to convert or reformat output files from three bioinformatics tools, Tandem Repeats Finder [(Benson, 1999)](https://www.ncbi.nlm.nih.gov/pubmed/9862982), WindowMasker from NCBI [(Morgulis et al., 2006)](https://www.ncbi.nlm.nih.gov/pubmed/16287941) and check_primers from Primer3 [(Rozen and Skaletsky, 2000)](https://www.ncbi.nlm.nih.gov/pubmed/10547847). All scripts have been tested with Chlamydomonas-specific datasets and may require some editing to become more widely applicable. 

# 1. BACSearcher
BACSearcher is a genome analysis python tool to facilitate easy recombineering of Chlamydomonas genes from BACs. 
Script: BACSearcher.py

### 1.1 Function
BACSearcher considers all or a subset of Chlamydomonas genes and provides three key pieces of information useful for recombineering a gene of interest using the pipeline detailed in Supplementary Method 1:
(1)	The script identifies up to five 50 bp regions chosen from a stretch of sequence 2000-3000 bp upstream of the annotated start codon for each gene. Regions are picked based on low GC content and absence of mono- and dinucleotide repeats. The script reports the reverse complement of each of these regions, any of which can be added to the 5’ end of the universal primer GAAGATCCTTTGATCTTTTCTACGGG to produce a target-specific cloning primer. When used in tandem with the 50 bp immediately upstream of the stop codon, prepended to the 5’ end of universal primer GGAGATCTGGGTGGCTCCG, these cloning primers can amplify a target-specific recombineering cassette from pLM099 by PCR (see Supplementary Method 1C). Instructions to modify the lengths of these regions are detailed below (see Modifications to the script). 
(2)	The script reports four pairs of Primer3-generated checking primers for each gene, two which can amplify a short sequence from the 5’ end of the gene and two from the 3’ end. These can be used to confirm the presence of the 5’ and 3’ ends of each target gene within a BAC in the Chlamydomonas BAC library.
(3)	The script reports the five smallest BACs containing each gene where available. BACs are reported if they cover the region spanning from 3000 bp upstream of the gene to the stop codon. An exception is made in cases where this region extends beyond the ends of a chromosome or scaffold, in which case the region is altered so that it starts/ends at the first/last position of the chromosome/scaffold. These parameters can be modified by the user so that BACs are only reported if, for example, the BAC sequence covers a larger 5’ flanking region, or also covers the 3’UTR of the gene of interest. Instructions to make these modifications to the script are detailed below. 

Two outputs are generated after each run of the script, one for BAC coverage and one for fosmid coverage of each gene. BAC and fosmid data are taken from the version 5.5 annotation of the Chlamydomonas genome. 
The fosmid output file contains identical information for (1) and (2), but for (3) the BAC-specific information is replaced by fosmid-specific information. Fosmid plate/well coordinates are not provided. The BAC output for all 17,741 genes in the genome is provided in Supplemental Data Set 1.

### 1.2 Specialist python modules required
BACSearcher has been tested using Python 3.6 and requires the following modules to be installed:
Gffutils	      https://pypi.org/project/gffutils/
Interval tree	  https://pypi.org/project/intervaltree/
Biopython	      https://pypi.org/project/biopython/
Primer3-py	    https://pypi.org/project/primer3-py/


### 1.3 Example usage
BACSearcher can be initiated from the command line using the following options:

`./BACSearcher.py 
-p BACs_fosmids.pairs.tsv 
-f Creinhardtii_281_v5.0.fa.gz 
-g Creinhardtii_281_v5.5.gene.gff3.gz 
-l Gene_shortlist.txt 
-w Bac_wells.txt 
-o Chlamydomonas_BACSearcher_results`

Where `BACs_fosmids.pairs.tsv` refers to precursor file I (below); `Creinhardtii_281_v5.0.fa.gz` refers to precursor file II; `Creinhardtii_281_v5.5.gene.gff3.gz` refers to precursor file III; `Gene_shortlist.txt` refers to precursor file IV; and `BAC_wells.txt` refers to precursor file V. `Chlamydomonas_BACSearcher_results` is the output file name. 

### 1.4 Required precursor files
* I:  TSV file containing the coordinates of the start and end of each valid BAC in the library. BACs are included as valid if their start and end sequences are mapped to the same chromosome and are in the correct orientation, (i.e., one end on each strand). This file is provided in the Supplemental Code ZIP folder as BACs_fosmids.pairs.tsv, is reproduced in Supplemental Data Set 5 and is also available from the GitHub repository (see Code availability section above).
* II:	Zipped FASTA file (.fa.gz extension) containing the gene sequences for all Chlamydomonas nuclear genes.
* III:	Zipped GFF file (.gff3.gz extension) containing version 5.5 annotation information for the Chlamydomonas genome.
* IV:	(Optional) TXT file containing the Cre IDs for all genes of interest to be processed, one per line, each appended with ‘.v5.5’. If this file is not provided, BACSearcher will process all nuclear genes and produce a TSV file of the results with the name specified by -o (see Example usage, above).
* V:	TXT file containing the plate and well coordinates of each BAC in the library, in the format ‘A-B-C’, where A is the plate number, B the row number and C the column number. This file is provided in the Supplemental Code ZIP folder as BAC_wells.txt and is also available from the GitHub repository (see Code availability section).
...VI.	(Optional) DB file generated from III using the BACSearcher script, which can be used in place of III in future runs. 

Genome FASTA and gene annotation GFF files (precursors II and III) are available for download from Phytozome, entitled Creinhardtii_281_v5.0.fa.gz and Creinhardtii_281_v5.5.gene.gff3.gz. The output provided in Supplemental Data Set 1 used precursor files II and III downloaded from Phytozome V12. 

When supplied with a GFF file via -g, BACSearcher will generate a gffutils database for the GFF file (precursor file VI). BACSearcher can also use this database directly with the -d option, saving the effort of regenerating the database.


### 1.5 Modifications to the script
The BACSearcher output for all 17,741 Chlamydomonas genes is provided in Supplemental Data Set 1 according to default parameters described in the Function section, above. Users can modify the script to change the length or position of the homology regions that BACSearcher reports, as well as the region that reported BACs should cover for each gene. These modifications are made by adding the additional options `-q`, `-r`, `-s`, `-t`, `-u` and `-v` to the command line before running the script.

1. Modifying the lengths of the homology arms:

...By default, BACSearcher reports 5’ and 3’ homology regions for each gene that are 50 bp long. To change these lengths to a different value, `x`, the options `-q` and `-r` can be added to the command line:

...`-q x` will change the default lengths of the reported 5’ homology regions to an integer, `x` bp.

...`-r x` will change the default length of the 3’ homology region to an integer, `x` bp.

...We recommend using the same values for `-q` and `-r`.

2. Modifying the size of the upstream native promoter region:

...By default, BACSearcher searches for 5’ homology regions 2000-3000 bp upstream of the start codon. To change the searched region, the default maximum flank of 3000 can be changed to a different value, `x`, by adding the option, `-s`, to the command line:

...`-s x` will direct the script to search for suitable homology regions in the 1000 bp downstream of an upstream position, `x` bp. For example, if `x=5000` the script will search for regions between 4000 and 5000 bp upstream of the start codon of each gene.

...If you would like to measure from the start of the 5’UTR instead of from the start codon, the option `-u` can be added to the command line:

...`-u` will direct the script to search for homology regions in the region defined by `-s` but measured from the 5’UTR instead of the start codon. If `-u` is used but `-s` is left undefined, the script will search for homology regions 2000-3000 bp upstream of the 5’UTR.

####	Modifying the BACs reported for each gene:
By default, the script is set to report BACs that cover the coding region plus the upstream flank for each gene, ignoring the 3'UTR. The upstream flank is defined by `-s` and `-u` (see (2), above) with a default value of 3000 bp upstream from the start codon. If you would like to add a downstream flank of length `x` bp that any reported BACs should also cover, the option `-t` can be added to the command line:
	`-t x` will direct the script to report only those BACs that cover `x` bp downstream from the stop codon of each gene.
If you would like to measure the downstream flank from the end of the 3’UTR instead of from the stop codon, the option `-v` can be added to the command line:
	`-v` will direct the script to report BACs that cover a downstream region defined by `-t` but measured from the 3’UTR instead of from the stop codon of each gene. If `-v` is used but `-t` is left undefined, the script will report BACs that cover up to the 3’ end of the 3’UTR.




Supplemental Code

Developed by Tom Emrich-Mills and John W. Davey

Associated manuscript: A recombineering pipeline to clone large and complex genes in Chlamydomonas, by Emrich-Mills et al. (2020)


BAC_wells.txt

Precursor file required by BACSearcher in order to provide plate coordinates for any predicted BACs within the Chlamydomonas BAC library.
For usage instructions, see Supplemental Method 2A; this file corresponds to precursor file V. 
The file is also available from the GitHub repository associated with this publication at the following link:


BACs_fosmids.pairs.tsv

Precursor file required by BACSearcher to predict BAC coverage for each gene.
For usage instructions, see Supplemental Method 2A; this file corresponds to precursor file I.
The file contains start and end coordinates for all valid BACs and fosmids within the Chlamydomonas BAC library, as well as the lengths of each. 
BACs are identified by a PTQ number and fosmids are identified by a VTP number. 
Plasmids are defined as valid if their start and end sequences are mapped to the same chromosome, with one end mapped to the sense strand and the other to the anti-sense strand. 
These criteria exclude 3179 of the 11,676 BACs (27.2%) and 20,696 of the 56,276 fosmids (36.8%) that are annotated against the Chlamydomonas genome (v5.5). 
Due to the way the BACs have been annotated against the genome, some valid plasmids have more than two ends mapped to the same chromosome and so are represented more than once (see Supplemental Method 2B for additional information).
This file is available from the GitHub repository associated with this publication at the following link:


BACSearcher.py

The BACSearcher python script, including brief usage instructions.
For more detailed information regarding function, usage and precautions, see Supplemental Method 2.
BACSearcher is also available from the GitHub repository associated with this publication at the following link:
https://github.com/TZEmrichMills/BACSearcher


Check_primers-output-sorter.txt

Python script intended to enable quick analysis of the output from the Primer3 check_primers module (Rozen and Skaletsky, 2000). 
Usage instructions are included at the top of the script and in Supplemental Method 4C. 
See Supplemental Method 3D for details of the Primer3 check_primers module.
This script is also available from the GitHub repository associated with this publication at the following link:


TRF-output-sorter.txt

Python script intended for processing the output information produced by Tandem Repeats Finder (Benson, 1999) into CSV format. 
Usage instructions are included at the top of the script and in Supplemental Method 4A.
See Supplemental Method 3A for details on Tandem Repeats Finder usage.
This script is also available from the GitHub repository associated with this publication at the following link:


WindowMasker-output-sorter.txt

Python script intended for counting the number of repeat sequences masked by WindowMasker (Morgulis et al., 2006).
The script takes the interval TXT file produced by WindowMasker and outputs a CSV file containing the number or repeats per input sequence.
Usage instructions are included at the top of the script and in Supplemental Method 4B.
See Supplemental Method 3C for details on WindowMasker and its usage.
This script is also available from the GitHub repository associated with this publication at the following link:

-
