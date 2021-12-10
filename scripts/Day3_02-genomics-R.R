# 02-genomics-R
# Day3 lesson 2

setwd("Bioinformatics_workshop-Dec-2021/Day-3/data")
##############################
library(org.Hs.eg.db)
org.Hs.eg.db
# what class is it
class(org.Hs.eg.db)
# what types of data can be extracted?
keytypes(org.Hs.eg.db)
# obtain all ENSEMBL IDs
entrez.ids <- keys(org.Hs.eg.db, keytype="ENSEMBL")
head(entrez.ids)

# how many are there
length(entrez.ids)
# use ensembl id of first 6 in entrez.ids to get desired keytypes
select(org.Hs.eg.db, keys = head(entrez.ids), columns = c("SYMBOL","ENTREZID", "REFSEQ"), keytype="ENSEMBL")
# using mapIds but only to get gene symbol
mapIds(org.Hs.eg.db, keys = head(entrez.ids), column = c("SYMBOL"), keytype="ENSEMBL")

# read in data
results <- read.csv("diff-exp-results.csv", stringsAsFactors = F, row.names = "ensembl")

# check the first few lines
head(results)

# using mapIds but only to get gene symbol
gene.symbols <- mapIds(org.Hs.eg.db, keys = rownames(results), column = c("SYMBOL"), keytype="ENSEMBL")

# add a symbols column to results with the gene symbols we just retreived
results$symbol <- gene.symbols

# check the new results table
head(results)

# make sure there are no NAs in the symbols column we just created
table(is.na(results$symbol))


library(EnsDb.Hsapiens.v86)


# using mapIds but only to get gene symbol
gene.symbols.2 <- mapIds(EnsDb.Hsapiens.v86, keys = entrez.ids, column = c("SYMBOL"), keytype="GENEID")

# how long is it
length(gene.symbols.2)

# how many NAs
table(is.na(gene.symbols.2))


anno <- read.table("GRCh38.p12_ensembl-97.txt", sep="\t", header=TRUE, stringsAsFactors = F)

# check the first few rows and dimensions
head(anno)
dim(anno)

# check how many Ensembl IDs overlap with our results
table(anno$Gene.stable.ID %in% rownames(results))
table(rownames(results) %in% anno$Gene.stable.ID)

# lets rename the ensembl ID column in both datasets so that we can merge them together based on those IDs
colnames(anno)[1] <- "ensembl"
results$ensembl <- rownames(results)

results_merge <- merge(results, anno, by="ensembl")
head(results_merge)
table(is.na(results_merge$Gene.name))

write.csv(results_merge, file = "diff-exp-results-annotated.csv")

#########################
# load the package
library(BSgenome)

# check currently available genomes
available.genomes()

# assign the genome to a variable using getBSgenome() (you need to have the package for the BSgenome you are trying to load already installed)
genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
genome

# check the structure
str(genome)

# print the metadata for the genome
metadata(genome)

genome.m <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")
class(genome.m)
genome.m

# unmasked genome
class(genome)
genome

# return basic sequence information summary
seqinfo(genome.m)

# print chromosome 1
genome.m$chr1

# assign chr 1
chr1 <- genome$chr1

# what is the frequency of each base in your sequence
alphabetFrequency(chr1, baseOnly=TRUE, as.prob=TRUE)

# what is the frequency of your favourite base
letterFrequency(chr1, "A", as.prob=TRUE)

# where are all incidences of 'ATG'
matchPattern("ATG", chr1)

# we need to establish a vector describing what the extra extended BED columns are
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

# read in peaks
bed <- import("CTCF-forebrain-mm10.bed",
              format="BED",
              extraCols = extraCols_narrowPeak,
              genome = "mm10")

# extract sequences for peak regions and print to console
ctcf_seqs <- getSeq(genome, bed)
ctcf_seqs

# calculate nucleotide freqs.
nt_freqs <- alphabetFrequency(ctcf_seqs, baseOnly=TRUE, as.prob=TRUE)

# calculate mean value for nucleotide freqs across all peaks
round(apply(nt_freqs, 2, mean), digits=2)

hist(width(ctcf_seqs),
     col = "darkgray",
     xlab = "Peak width (bp)",
     main = "CTCF peak width distribution")

# resize the regions from the BED file
bed_centered <- resize(bed, width = 400, fix = "center")
bed_centered

# check their with
width(bed_centered)

# extract sequences again
ctcf_seqs_cent <- getSeq(genome, bed_centered)
ctcf_seqs_cent

# add names to peaks in ctcf_seqs so that FASTA entries have names
names(ctcf_seqs) <- paste0(seqnames(bed), ":", start(bed), "-", end(bed))

# export peaks to FASTA file
writeXStringSet(ctcf_seqs, file="CTCF-peaks-resized.fa")
