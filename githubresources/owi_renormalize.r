#! /usr/bin/env Rscript

## Script for renormalizing abundances of a csv file
## The script will read the following arguments from the command line: the input file name, the output file name,
## the threshold (cutoff) which will be used to equal to zero all cumulative frequencies below that cutoff
## and the numbers for the starting and end sample columns 
## If the output file argument is empty, it will just add ".renormalized.csv" at the end of the name of the input file.
## By Owen S. Wangensteen - Project Metabarpark  2015

library("optparse")

option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="Input tsv file with abundance columns", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL, 
              help="Output file name [default = input file ending in renormalized.tsv]", metavar="character"),
  make_option(c("-c", "--cutoff"), type="numeric", default=0.99, 
              help="Threshold for trimming [default = 0.99]", metavar="numeric"),
  make_option(c("-s", "--start"), type="integer", default= 14, 
              help="Sample columns start [default = 14]", metavar="numeric"),
  make_option(c("-e", "--end"), type="integer", default= 88, 
              help="Sample columns end [default = 88]", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$infile) ){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input tsv file with abundances).", call.=FALSE)
}

if (is.null(opt$outfile)){
  opt$outfile <- paste(substr(opt$infile,1,nchar(opt$infile)-3),"renormalized.tsv",sep="")
}
print(opt)
showlines <- 500
min_abund <- 99 # Minimum total abundance for a MOTU to be renormalized

# Read Database
  message("Reading database...")
  db <- read.table(opt$infile,sep="\t",head=T,quote = "\"") #only look for double quotes!
  if (ncol(db)==1) db <- read.table(opt$infile,sep=",",head=T) # replace by ; if , is used for decimals, dec=
  total_seqs <- nrow(db)
  sample_columns <- opt$start:opt$end
  all_numeric = 0
  for (i in sample_columns)  all_numeric <- all_numeric + (class(db[,i])=="integer")
  if (all_numeric != length(sample_columns)) stop("Some columns defined as samples include non integer data! Please check start and end options.")
  initial_reads <- sum(db[,sample_columns])
  initial_counts <- rowSums(db[,sample_columns])
  message("Database read including ", total_seqs," total MOTUs and ",initial_reads," initial reads.")
  
  #Extract counts database
  counts.db <- db[,sample_columns]
  
  #Calculate frequency matrix
  message("Calculating renormalized matrix.")
  suma.samples <- colSums(counts.db)
  freq.matrix <- counts.db #For dimensioning
  for (columna in 1:ncol(counts.db)) freq.matrix[,columna] <- counts.db[,columna]/suma.samples[columna]
  check <- colSums(freq.matrix)
  
  #Select MOTUs to renormalize
  suma.motus <- rowSums(counts.db)
  to_renormalize <- suma.motus > min_abund
  message(sum(to_renormalize)," MOTUs with abundance > ",min_abund," will be renormalized.")
  
  #Calculate renormalized freq matrix
  
  suma.taxa <- rowSums(freq.matrix)
  freq.matrix.normalized <- freq.matrix #For dimensioning
  for (fila in (1:nrow(freq.matrix))[to_renormalize]) freq.matrix.normalized[fila,] <- freq.matrix[fila,]/suma.taxa[fila]
  check2 <- rowSums(freq.matrix.normalized)
  
  #Delete all counts below cutoff of cumulative freq
  counts.trimmed <- counts.db
  for (fila in (1:nrow(freq.matrix))[to_renormalize]){
    numeros <- as.vector(freq.matrix.normalized[fila,order(freq.matrix.normalized[fila,],decreasing=T)],"numeric")
    acumulado <- cumsum(numeros)
    muestra.limit <- sum(acumulado<opt$cutoff)+1
    threshold <- numeros[muestra.limit]
    counts.trimmed[fila,][freq.matrix.normalized[fila,]<threshold] <- 0
    if (i %% showlines == 0) message(fila,"/",nrow(freq.matrix)," MOTUs calculated.","\r",appendLF = FALSE)
  }
  
  #Reconstruct the database
  #new.db <- data.frame(db[,1:opt$start-1],counts.trimmed,initial_counts,final_counts=rowSums(counts.trimmed),counts_trimmed=initial_counts-rowSums(counts.trimmed),db[,(opt$end+1):ncol(db)])
  new.db <- data.frame(db[,1:opt$start-1],counts.trimmed,initial_counts,final_counts=rowSums(counts.trimmed),counts_trimmed=initial_counts-rowSums(counts.trimmed))
  new.db_escasos <- new.db[new.db$final_counts<=2,]
  new.db_singletons <- new.db[new.db$final_counts<=1,]
  
  #Save the final database
  write.table(new.db,opt$outfile,sep="\t",row.names=F)
  final_reads <- sum(new.db$final_counts)
  message("File ", opt$outfile, " written, including ",final_reads," reads in ",nrow(new.db)," MOTUs.")
  message(sum(new.db$counts_trimmed)," total reads deleted while trimming.")
  

