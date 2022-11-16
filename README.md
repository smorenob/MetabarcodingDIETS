# MetabarcodingDIETS
 
This workflow is specific to analysing eukaryote diversity from community DNA metabarcoding of gut contents from different species at zoantharian dominated habitats using CO1 gene region from raw reads. It reproduces the cleaning and curation steps as well as the biostatiscal analysis included in the manuscript "Metabarcoding hyperdiverse kelp holdfast communities on temperate reefs: an experimental approach to inform future studies." by Vanessa Arranz, Libby Liggins and J. David Aguirre.  

Sample preparation, DNA extraction and PCR amplification steps related to this work can be found here (add link).

The scripts are designed to be run using a Linux OS, and were developed on Ubuntu 16.04. 
#insert line

## Requirements : 

- The contents of this repository
- Raw sequence data : All FASTQ sequence files are available from the National Center for Biotechnology Information short-read archive database (Bioproject: PRJNA638997, Biosamples: SAMN15220525-SAMN15220620).
- QIIME2 version 2022.8 https://docs.qiime2.org/2020.8/install/ 
- Biom http://biom-format.org/
- MARES reference sequences database : https://osf.io/4f8mk/ 
- BLASTn  https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
- MEGAN6 -  Metagenome Analyzer https://www.wsi.uni-tuebingen.de/lehrstuehle/algorithms-in-bioinformatics/software/megan6/
- VSEARCH https://github.com/torognes/vsearch
- R https://www.r-project.org/

## Getting started

### Activate QIIME2 

**Citation:** Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). https://doi.org/10.1038/s41587-019-0209-9

```
conda activate qiime2-2022.8
``` 
### Import raw sequence data


It is important to check what input-format we are going to use. In our case it is a FastQ with one file (forward and reverse are analyzed separately in this case). Files are splitted into one different folder for each sample an inside there will be the file. If that it is the case, we need to extract all files into a single folder before importing. 

For forward files:
```

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path zoandietFfastq/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-single-endF.qza
```  
For reverse files:
```
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path zoandietRfastq/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-single-endR.qza

# paired-end-demux.qza is an artifact file '.qza' that has the information stored for your raw fastq sequences. We can summarize and visualize the output by transforming to '.qzv'. Check the output : https://view.qiime2.org/

qiime demux summarize \
  --i-data demux-single-endF.qza \
  --o-visualization demux-single-endF.qzv

qiime tools view demux-single-endF.qzv
```

## Remove primers : cutadapt

We use [cutadapt](https://github.com/qiime2/q2-cutadapt) to remove the primers

**Citation:** Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 17(1):pp–10, 2011. https://doi:10.14806/ej.17.1.200.

```
qiime cutadapt trim-single \
  --i-demultiplexed-sequences demux-single-endF.qza \
  --p-front GGWACWGGWTGAACWGTWTAYCCYCC \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-seqsF.qza \
  --verbose
  
qiime demux summarize \
  --i-data trimmed-seqsF.qza \
  --o-visualization trimmed-seqsF.qzv
  
qiime tools view trimmed-seqsR.qzv

qiime cutadapt trim-single \
  --i-demultiplexed-sequences demux-single-endR.qza \
  --p-front TANACYTCNGGRTGNCCRAARAAYCA \
  --p-error-rate 0 \
  --o-trimmed-sequences trimmed-seqsR.qza \
  --verbose
  
qiime demux summarize \
  --i-data trimmed-seqsR.qza \
  --o-visualization trimmed-seqsR.qzv
  
qiime tools view trimmed-seqsR.qzv
```
## Denoise, chimera removal and clustering into ASVs: DADA2

[DADA2](https://github.com/qiime2/q2-dada2) : Pair-end joining, dereplication, chimera filtering and clustering in ASVs 

**Citation:** Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869

We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments). (In this case this does not applie since we are working with Forward and Reverse separately)

We can now construct an Amplicon Sequence Variant(ASV) table.

ASVs are inferred by a *de novo* process in which biological sequences are discriminated from errors on the basis of, in part, the expectation that biological sequences are more likely to be repeatedly observed than are error-containing sequences. As a result, ASV inference cannot be performed independently on each read—the smallest unit of data from which ASVs can be inferred is a sample. However, unlike *de novo* OTUs, ASVs are consistent labels because ASVs represent a biological unit that exists outside of the data being analyzed. Thus, ASVs inferred independently from different studies or different samples can be validly compared.

```
qiime dada2 denoise-single \
--i-demultiplexed-seqs trimmed-seqsF.qza \
--p-trunc-q 30 \
--p-trunc-len 125 \
--p-chimera-method consensus \
--o-table table-dada2F.qza \
--o-representative-sequences rep-seqs-dada2F.qza \
--o-denoising-stats denoising-stats-dada2F.qza 

qiime feature-table summarize \
--i-table table-dada2F.qza \
--o-visualization table-dada2F.qzv 

qiime tools view table-dada2F.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2F.qza \
--o-visualization rep-seqs-dada2F.qzv 

qiime tools view rep-seqs-dada2F.qzv

# From this visualization we can download the FASTA FILE that we will use for taxonomic assignment (see Taxonomic assignment section)

qiime metadata tabulate \
--m-input-file denoising-stats-dada2F.qza \
--o-visualization denoising-stats-dada2F.qzv

qiime tools view denoising-stats-dada2F.qzv
```
## Export ASV table and representative sequences
```
# Example ASV table without filtering
qiime tools export \
--input-path table-dada2F.qza \
--output-path ASVFtable/

#Export representative sequences
qiime tools export  \
--input-path rep-seqs-dada2F.qza \
--output-path rep-seq-ASVF.fasta
```
Convert the ASV table to tsv in Biom
```
biom convert \
-i ASVFtable/feature-table.biom \
-o ASVFtable/ASV-frequency-table.tsv  --to-tsv
```
##Taxonomic assignment
ASVs passing the quality control and filtering thresholds (rep-seq-ASV.fasta) were taxonomically assigned using the MARES reference sequence database. (In this case it was MIDORI2)

MARES is the most comprehensive CO1 reference database for marine eukaryotes available, and provides users the ability to retain taxa that cannot be assigned at the species level, but can be assigned at higher taxonomic levels.

To use MARES reference database:

Download MARES_NOBAR from https://osf.io/4f8mk/ Note : You can also use MARES_BAR.

Place MARES_NOBAR_BOLD_NCBI_sl_reformatted.fasta in the main folder or include the file path after -db

Build the database
We first create a blast database from the MARES reference sequence database
```
makeblastdb -in MARES_NOBAR_BOLD_NCBI_sl_reformatted.fasta -dbtype nucl -parse_seqids
```
Blast each sequence against the MARES reference sequences database : Blastn
Then, we performed a BLASTn against MARES reference database with an e-value of 1-60 for high-quality matches and with default max_target_seqs (500).
```
blastn -db midori2/MIDORI2_UNIQ_NUC_GB251_CO1_BLAST.fasta -query rep-seq-ASVF.fasta/dna-sequences.fasta -evalue 1e-60 -outfmt 5 -out MIDORIF_MEGAN.txt -num_threads 8
```
##Use LCA algorithm for taxonomic assignment : MEGAN6

We used MEGAN 6.18.9 for taxonomic assignment within the NCBI taxonomy framework using the default Lowest Common Ancestor (LCA) algorithm parameters.

**Citation** : Huson DH, Beier S, Flade I, Górska A, El-Hadidi M, et al. (2016) MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data. PLOS Computational Biology 12(6): e1004957. https://doi.org/10.1371/journal.pcbi.1004957

Launch MEGAN 6.24.1

File -> Import: the blast output and the fasta file used as the blast query into MEGAN (File → Import From BLAST) : MARES_MEGAN.txt and rep-seq-ASV.fasta

Apply the following LCA settings (Options -> LCA settings):
```
min score 98 
max expected 0.00000001 
min % ID 75 
top % 10 
min support % 0 (off) 
min support 1 
```
Select level of taxonomy to view, e.g. species, genus, family (Options -> Project assignment...). 

File -> Export -> Text (csv) Format Choose: readName_to_taxonPathKPCOFGS

- The .csv file will have a percent value at the end of the line. This refers to the percentage of high scoring alignments for the given read that map to the last taxon on the path. It has nothing to do with the percentage used in the weighted LCA.
- It only reports taxa in the path that have an official KPCOFGS rank. Intermediate nodes that have no taxonomic rank, or one that does not belong to KPCOFGS, are suppressed KPCOFGS = Kingdom, Phylum, Class, Order, Family, Genus, Species
- Each node is prefixed by a letter__ to indicate the rank, e.g. g__ for genus, s__ for species : we can edit this later to have only the name

Save into the main folder as assigned_seqs-MARES-ex.txt

**To create the OTU_taxonomy file to use in downstream analysis (see Statistical analysis):

Export it as ReadName_to_taxonPathKPCOFGS / ASSIGNED
Import it into excel and remove the columns with the percentage *No need in MEGAN6.21.10
Because the unknowns create problems, the excel spreadsheet needs to be edited to replace them with the taxon known
The ASVs not assigned have to be set as 'd_unassigned' because is not working otherwise
Generate the names with a semicolon ; separation
Create OTU_ID in the first column and the TaxonPath in the other one
Save as .csv or .txt -> taxonomy_edited_8ranks.txt
