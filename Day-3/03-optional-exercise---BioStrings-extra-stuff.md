# Optional Exercise: Additional *BioStrings* functionality

This optional exercise accompanies lesson *Working with genomics data in R/Bioconductor - Part III (Sequence analysis & references genomes)*. You may complete this exercise at home after the workshop, or during any session where you have extra time.

---

### Using *BioStrings* without *BSGenome*

BioStrings can be used independently from BSGenome with any set of sequences you are able to define in your R environment as an *XString* or *XStringSet* class object. For example, perhaps you are studying the *Amphimedon queenslandica*, a marine sponge organism native to the Great Barrier Reef, and want to explore some basic features of its coding sequences.

<p align="center">
<img src="../figures/Amphimedon_queenslandica_adult.png" title="xxxx" alt="context"
	width="70%" height="70%" />
</p>

Image source: [Wikipedia](https://en.wikipedia.org/wiki/Amphimedon_queenslandica)

We can retrieve a FASTA file for the coding sequences  (13 overall) from NCBI [(RefSeq ID: NC_008944.1)](https://www.ncbi.nlm.nih.gov/genome/2698) and read the FASTA file into R as a DNAStringSet object using the `readDNAStringSet()` function.
```{r}
fasta.file <- "a.queenslandica.fasta"
a.queen <- readDNAStringSet(fasta.file, "fasta")
a.queen
```

Just as we have done earlier in this lesson, we can again use the BioStrings functions to perform basic operations on these sequences. For example:
```r
# confirm how long it is
length(a.queen)

# what is the frequency of each base in your sequence
base.freqs <- alphabetFrequency(a.queen, baseOnly=TRUE, as.prob=TRUE)
base.freqs

# what is the frequency of your favourite base
a.freqs <- letterFrequency(a.queen, "A", as.prob=TRUE)
a.freqs
```

BioStrings also implements extensive functionality for **pattern matching**, allowing you to search sequences for specific patterns of interest. For example, we may want to confirm that each of our coding sequences begins with an `ATG` start codon. We can do this using the BioStrings functions `matchPattern()` and `countPattern()`.
```r
# return all matches in a DNAString subject sequence to a query pattern
matchPattern("ATG", a.queen[[1]])

# only return how many counts were found
countPattern("ATG", a.queen[[1]])

# what happens if we remove the indexing of the DNAStringSet object 'a.queen'? Why?
matchPattern("ATG", a.queen)

# match a query sequence against multiple sequences that you want to search in  
vmatch <- vmatchPattern("ATG", a.queen)
vmatch

# look at the structure of this object
str(vmatch)

# extract the IRanges for matches in the first subject sequence
vmatch[[1]]
```

We may also have several patterns that we want to search for in each our coding sequences. For example, perhaps we want to search for standard stop codons (`TAG`, `TAA`, `TGA`) in the *A.queenslandica* coding sequences. BioStrings functions `matchPDict()` and `vmatchPDict()` provide functionality for such tasks. e.g.
```r
# create a DNAStringSet of the stop codons
stop.codons <- DNAStringSet(c("TAG", "TAA", "TGA"))

# create a dictionary of patterns (PDict class object) that you want to search your sequences for
stop.codons.dict <- PDict(stop.codons)

# search in the first coding sequence
match1 <- matchPDict(stop.codons.dict, a.queen[[1]])
match1
match1[[3]]

# use a loop to search for stop codons in all our coding sequences
matches <- list()
for(i in 1:length(a.queen)){
  matches[[i]] <- matchPDict(stop.codons.dict, a.queen[[i]])
}
length(matches)
str(matches)
matches[[4]]
matches[[4]][[3]]
```

---

#### Other functionality in BioStrings

BioStrings also provides functionality for a number of other analytical tasks that you may want to perform on a set of sequences stored using the *XString* and *XStringSet* method, for example:  
* trimming sequence ends based on pattern matching using `trimLRPatterns()`
* local and global alignment problems using `pairwiseAlignment()`
* read in multiple sequence alignments using `readDNAMultipleAlignment()`
* motif searches with a Position Weight Matrix (PWM) using `matchPWM()` (commonly done in ChIP-seq & ATAC-seq)
* palindrome searching using findPalindromes `findPalindromes()`
* computing edit distances between sets of sequences using `stringDist()`  
