# MetabarcodingDIETS
 
This workflow is specific to analysing eukaryote diversity from community DNA metabarcoding of gut contents from different species at zoantharian dominated habitats using CO1 gene region from raw reads. It reproduces the cleaning and curation steps as well as the biostatiscal analysis included in the manuscript "Metabarcoding hyperdiverse kelp holdfast communities on temperate reefs: an experimental approach to inform future studies." by Vanessa Arranz, Libby Liggins and J. David Aguirre.  

Sample preparation, DNA extraction and PCR amplification steps related to this work can be found here (add link).

The scripts are designed to be run using a Linux OS, and were developed on Ubuntu 16.04. 
#insert line

## Requirements : 

- The contents of this repository
- Raw sequence data : All FASTQ sequence files are available from the National Center for Biotechnology Information short-read archive database (Bioproject: PRJNA638997, Biosamples: SAMN15220525-SAMN15220620).
- QIIME2 version 2019.4 https://docs.qiime2.org/2020.8/install/ 
- Biom http://biom-format.org/
- MARES reference sequences database : https://osf.io/4f8mk/ 
- BLASTn  https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- MEGAN6 -  Metagenome Analyzer https://www.wsi.uni-tuebingen.de/lehrstuehle/algorithms-in-bioinformatics/software/megan6/
- VSEARCH https://github.com/torognes/vsearch
- R https://www.r-project.org/
