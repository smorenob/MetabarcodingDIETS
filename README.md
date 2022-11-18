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

File -> Import: the blast output and the fasta file used as the blast query into MEGAN (**File → Import From BLAST**) : MARES_MEGAN.txt and rep-seq-ASV.fasta

Apply the following LCA settings (**Options -> LCA settings**):
```
min score 98 
max expected 0.00000001 
min % ID 75 
top % 10 
min support % 0 (off) 
min support 1 
```
Select level of taxonomy to view, e.g. species, genus, family (**Options -> Project assignment...**). 

**File -> Export -> Text (csv) Format Choose: readName_to_taxonPathKPCOFGS**

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

##Create ASV/OTU Phyloseq objects
We then import the ASV table in R as phyloseq object with the taxonomy associated, sample metadata and reference sequences for statistical analysis.
```
library("phyloseq")
library("readr")
library("dplyr")
library("Biostrings")

# ASV table without filtering (change format from .tsv to .csv)

ASV_table <- read.csv("ASVFtable/ASV-frequency-table.csv", row.names = 1, header = TRUE, check.names = FALSE)
ASV_table <- otu_table(as.matrix(ASV_table), taxa_are_rows = TRUE)

## Add the taxonomy assigned 

tax_table_new_edited_8ranks <- read.csv("taxonomy_edited_8ranksF.csv")
str(tax_table_new_edited_8ranks)

matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[,1]
  m
}
newtaxonomy <-matrix.please(tax_table_new_edited_8ranks)
tax_table_newtaxonomy_8ranks <- tax_table(newtaxonomy)

# Add the sample data file 
sample_data_96samples <- read.csv("Resources/sample_data_5X96samples.csv")
sampledata = sample_data(data.frame(sample_data_96samples, row.names = sample_names(ASV_table)))

# Add the reference sequences

reference_seqs0 <- readDNAStringSet(file = "rep-seq-ASVF.fasta/dna-sequences.fasta",format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

# Add the sample data file 
sample_data_96samples <- read.csv("Resources/sample_data_5X96samples.csv")
sampledata = sample_data(data.frame(sample_data_96samples, row.names = sample_names(ASV_table)))

# Create the phyloseq object 
ASV <- phyloseq(otu_table(ASV_table), sample_data(sampledata), refseq(reference_seqs0), tax_table(tax_table_newtaxonomy_8ranks))
```
##Extract possible contaminants and tag switching normalisation
We used decontam package in R to identify and extract contaminants using the DNA concentration of the samples and the negative controls.

Citation: Davis NM, Proctor D, Holmes SP, Relman DA, Callahan BJ (2017). “Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data.” bioRxiv, 221499. https://doi.org/10.1101/221499
```
library("phyloseq")
library("decontam")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

# Identify and extract blank contaminants using the negative controls. The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants ###

sample_data(ASV)$is.neg <- sample_data(ASV)$Sample_or_Control == "Control"
contam.prev <- isContaminant(ASV, method="prevalence", neg="is.neg",  threshold=0.5)
table(contam.prev$contaminant)

# Extract contaminants
ASV.nocon <- prune_taxa(!contam.prev$contaminant, ASV)

# Extract negative control samples 
ASV.nocon.nc <- subset_samples(ASV.nocon, Sample_or_Control=="True Sample")

write.csv(otu_table(ASV_nf.noSC.nocon.nc),file = "ASVtable_nf_nosc_noc_nocon.csv")
```
Tag switching correction. We used the R Script included in Resources/owi_renormalize.R(https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast/blob/main/Resources/owi_renormalise.R). This sorts the samples by abundance for each ASV and eliminates the reads of the samples corresponding to a cumulative frequency of less than 3% for each particular ASV.

Citation: Wangensteen OS, Turon X (2016) Metabarcoding techniques for assessing biodiversity of marine animal forests. Marine Animal Forests. The Ecology of Benthic Biodiversity Hotspots, eds Rossi S, Bramanti L, Gori A, Orejas C (Springer International Publishing).
```
# In terminal (install package "optparse" in R):

RScript owi_renormalize.R -i ASVtable_nocon.csv -o ASVtable_tsc.csv -c 0.97 -s 2 -e 321
```
Rename samples names from . to - in the output table ASVtable_tsc.tsv. In R :
```
#After tag switching normalization
ASVtable_tsc_all <- read.csv("Resources/ASVtable_tsc.csv",sep = "\t", header=TRUE,as.is=TRUE, row.names = 1, check.names = FALSE)
str(ASVtable_tsc_all)
# delete the last 3 columns of total counts 
ASVtable_tsc_all <- ASVtable_tsc_all[-c(321:323)]
str(ASVtable_tsc_all)
# Remove rows that sum columns are 0
ASVtable_tsc_all <- ASVtable_tsc_all[as.logical(rowSums(ASVtable_tsc_all != 0)), ]
str(ASVtable_tsc_all)

ASV_nofilter <- otu_table(ASVtable_tsc_all, taxa_are_rows = TRUE)

ASV_nofilter <- phyloseq(otu_table(ASV_nofilter), sample_data(ASV.nocon.nc), refseq(ASV.nocon.nc), tax_table(ASV.nocon.nc))
```

##Clustering ASVs into OTUs : VSEARCH
We used Vsearch to cluster the ASVs into Operational Taxonomic Units (OTUs). We chose 97% of similarity for CO1 mitochondrial region.

Citation: Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584. https://doi.org/10.7717/peerj.2584
```
vsearch --cluster_size rep-sequences-nofilt.fasta  --id 0.97 --uc clustering-results.uc -msaout sequences-outs

qiime vsearch cluster-features-de-novo \
  --i-table table.qza \
  --i-sequences rep-seqs.qza \
  --p-perc-identity 0.99 \
  --o-clustered-table table-dn-99.qza \
  --o-clustered-sequences rep-seqs-dn-99.qza
```
Save the clustering results.uc as .csv

Manually deleting the Cluster records rows "C", because they are the same as the cluster centroids "S" -> clustering-results.csv

Change the header names to ASV_ID to merge in R with the previous taxonomy assigned to ASVs.

By running the following script in R, we created the OTU table by combining the previous ASV table and Vsearch clustering results using the ASV_ID as an indicator.
```
clustering.results <- read.csv("Resources/clustering-results.csv")
ASV_OTU_table <- merge(ASVtable_tsc_all_clean, clustering.results, by.x="ASV_ID", by.y="ASV_ID")
write.csv(ASV_OTU_table, file = "ASV_OTU_table.csv")
```
Manually edited again the ASV_OTU_table.csv to create the ASV_OTU_tabletocollapse.csv

Sort by Record type column : copy all the ASV_ID column of the centroids "S" into the OTU column (they are the consensus of the clusters)
Sort by Cluster number
Delete all the columns from Vsearch output + taxonomy and ASV_ID columns
```
OTU_nofilter_tocollapse <-  read.csv("Resources/ASV_OTU_tabletocollapse.csv", check.names = FALSE)
OTU_nofilter_tocollapse$OUT97_ID <- as.factor(OTU_nofilter_tocollapse$OUT97_ID)
#Collapse the OTUs with the same name (all ASVs that were collapse into OTUs)
OTU_nofilter_collapsed <- OTU_nofilter_tocollapse %>% group_by(OUT97_ID) %>% summarise_all(funs(sum))

OTU_nofilter_collapsed$OUT97_ID <- as.character(OTU_nofilter_collapsed$OUT97_ID)
OTU_nofilter_collapsed <- as.data.frame(OTU_nofilter_collapsed)
rownames(OTU_nofilter_collapsed) <- OTU_nofilter_collapsed[,1]
OTUtable_tsc <- OTU_nofilter_collapsed[,-1]


# Create Phyloseq object with OTU table after tag switching normalisation 

OTU_nofilter <- otu_table(OTUtable_tsc, taxa_are_rows = TRUE)
OTU_nofilter <- phyloseq(otu_table(OTU_nofilter), sample_data(sampledata00), refseq(reference_seqs0), tax_table(tax_table_newtaxonomy_8ranks))
```
