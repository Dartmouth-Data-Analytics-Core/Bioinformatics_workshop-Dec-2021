# Optional Exercise: Programmatic access of BioMart

This optional exercise accompanies lesson *Working with genomics data in R/Bioconductor - Part II (Genome Annotation)*. You may complete this exercise at home after the workshop, or during any session where you have extra time.

---

### Access BioMart data from within R

BioMart can be accessed programmatically from within R, using the `BioMart` package. This can be a slower than the approach described above, but allows more flexibility depending on the annotation data required.

```r
library(biomaRt)

# check available martslistMarts(), 10)
head(listMarts())
```

You can see the same marts are listed as were available from the BiomaRt website. We need to choose one of these marts in order to see the available annotation datasets that can be accessed from that mart.
```r
# use 'ensembl to select the 'ENSEMBL_MART_ENSEMBL' mart
ensembl <- useMart("ensembl")

# show available datasets to pull annotation data from
head(listDatasets(ensembl), 10)
tail(listDatasets(ensembl), 10)
nrow(listDatasets(ensembl))

# check for human dataset
table(listDatasets(ensembl)$dataset %in% "hsapiens_gene_ensembl")
```

Now select the `hsapiens_gene_ensembl` dataset from the mart and view the available *attributes* (values that can be returned, e.g. gene symbol) for that dataset.
```r

# pick the ensembl mart for humans  
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# list the attributes for this dataset
head(listAttributes(ensembl), 10)
tail(listAttributes(ensembl), 10)
nrow(listAttributes(ensembl))
```

The flagship function of the BiomaRt package is `getBM()` (for get BiomaRt presumably) which allows us to obtain specific data (attributes) given a set of values that we provide (e.g. Ensembl IDs). The process is very similar to how we used the `select()` and `mapIDs()` functions from OrgDb. Lets use `getBM()` to return annotation data from BiomaRt for our RNA-seq data.
```r
# submit query for 1st 2000 genes (to save time in class) in our RNAseq results
anno_bm <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
                 filters = "ensembl_gene_id",
                 values = head(results$ensembl, 2000),
                 mart = ensembl,
		 useCache = FALSE)
head(anno_bm, 10)
```

You can see we now have a `data.frame` stored in anno_bm in our environment that looks similar to the text file that we downloaded directly from biomaRt and read into R. You could follow a similar protocol to that which we performed to merge the BiomaRt data downloaded in the text file with our RNA-seq results.
