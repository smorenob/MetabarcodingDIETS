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
  
For reverse files:

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path zoandietRfastq/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-single-endR.qza
  ```

# paired-end-demux.qza is an artifact file '.qza' that has the information stored for your raw fastq sequences. We can summarize and visualize the output by transforming to '.qzv'. Check the output : https://view.qiime2.org/

```
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
  --p-front TANACYTCNGGRTGNCCRAARAAYCA \
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

We now merge the forward and reverse reads together to obtain the full denoised sequences. Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region (but these conditions can be changed via function arguments).

We can now construct an Amplicon Sequence Variant(ASV) table.

ASVs are inferred by a *de novo* process in which biological sequences are discriminated from errors on the basis of, in part, the expectation that biological sequences are more likely to be repeatedly observed than are error-containing sequences. As a result, ASV inference cannot be performed independently on each read—the smallest unit of data from which ASVs can be inferred is a sample. However, unlike *de novo* OTUs, ASVs are consistent labels because ASVs represent a biological unit that exists outside of the data being analyzed. Thus, ASVs inferred independently from different studies or different samples can be validly compared.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 200 \
  --p-trunc-len-r 200 \
  --p-n-threads 4 \
  –-p-chimera-method consensus \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza
  
qiime feature-table summarize \
--i-table table-dada2.qza \
--o-visualization table-dada2.qzv 
qiime tools view table-dada2.qzv

qiime feature-table tabulate-seqs \
--i-data rep-seqs-dada2.qza \
--o-visualization rep-seqs-dada2.qzv \
qiime tools view rep-seqs-dada2.qza

# From this visualization we can download the FASTA FILE that we will use for taxonomic assignment (see Taxonomic assignment section)

qiime metadata tabulate \
--m-input-file denoising-stats-dada2.qza \
--o-visualization denoising-stats-dada2.qzv
```
