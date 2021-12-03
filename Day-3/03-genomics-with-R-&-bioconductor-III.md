
# Working with genomics data in R/Bioconductor - Part III

## Sequence analysis & references genomes

Beyond providing access to annotation data in R, Bioconductor also provides functionality for accessing and analyzing complete reference sequences for commonly used genomes. Namely, the [BSgenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) family of Bioconductor packages provides an efficient way to obtain, query, and manipulate genomic sequence data from reference genomes. You can return a vector of the currently available genomes to your console by printing `available.genomes()` after loading the `BSgenome` package.

Analyzing genomic sequence data directly can be used for a number of common research tasks, all possible in the BSGenome/BioStrings framework, for example:  
* Extracting DNA/RNA/protein sequences for specific genomic features  
* Calculating nucleotide frequencies for defined sequences
* Searching for matching sequences of interest

```r
# load the package
library(BSgenome)

# check currently available genomes
available.genomes()
```

These genomes are focused predominantly on those available from NCBI and UCSC genomes, however functionality exists to [forge a BSGenome](https://bioconductor.org/packages/release/bioc/html/BSgenome.html) package, allowing you to leverage the BSGenome framework for genomes not part of the currently available set.

BSgenome packages are heavily dependent on the [BioStrings package](http://bioconductor.org/packÂages/release/bioc/html/Biostrings.html), which defines a set of methods and object classes for storing and accessing sequence data. BioStrings is loaded automatically when you loaded BSgenome.

---

## Learning objectives:

We will introduce the basic object classes and methods introduced by BioStrings to demonstrate how they form the basis for BSgenome packages, and what genomic sequence-based operations they allow you to perform on reference genome data (or any sequence you define).

---

### Basic object-classes and methods in the BioStrings package

Before working with a complete genome sequence from BSGenome, lets discuss the basic object-classes implemented in the BioStrings package to store sequence data, as well as the methods use to parse these sequences.

The most basic object class in BioStrings is the *XString* class, which is technically a '*virtual class*' (meaning it cannot actually store objects itself, but can be used to set rules for a a group of classes) encompassing object classes for DNA, RNA, and protein sequences: `DNAString`, `RNAString`, and `AAString`. Today we will focus on DNA sequences using DNAString class objects.

Lets start by creating a simple DNAString object and looking at some of its basic features:
```r
# use the DNAString constructor function to create a 10 letter DNA sequence
seq <- DNAString(x="AGCT", start=1, nchar=NA)
seq

# how long is it
length(seq)

# show the structure of the DNAString object
str(seq)
```

Now we will make a longer sequence using the pre-stored BioStrings object `DNA_ALPHABET` and some functions from base R.
```r
# print DNA alphabet to see what it returns
DNA_ALPHABET
```

The four standard DNA bases are returned as the first 4 elements of this character string. The remaining elements represent ambiguous bases or specific combinations/relationships using something called the *extended The International Union of Pure and Applied Chemistry (IUPAC) genetic alphabet.* BioStrings and its object classes use the extended IUPAC genetic alphabet to describe nucleic acid sequences, therefore we will briefly cover the basics of the extended IUPAC alphabet now.

The extended IUPAC code is a nomenclature used to describe incompletely specified nucleic acids. This is useful when dealing with reference genomes, as we often encounter regions with ambiguous assemblies. The standard IUPAC code uses 16 characters to specify both single bases (A, G, C, T, U) or the possible states for that base.

The standard IUPAC code is used by numerous bioinformatics tools and softwares in order to represent complex sequences of nucleic acids for which we may not be confident in some individual base identities (e.g. complex genomic regions that are challenging to sequence using short read approaches).

**Table 1. Standard IUPAC genetic alphabet.**

|Symbol |	Mnemonic| Translation  |
|---|---|---|
| A	|	A | (adenine) |  
| C	|	C | (cytosine)  |
| G	|	G	| (guanine)  |
| T	|		T	| (thymine)  |
| U	|		U	| (uracil)  |
| R	|	pu**R**ine		| A or G  |
| Y		| p**Y**rimidine		| C or T/U  |
| S		| **S**trong interaction	|	C or G  |
| W		| **W**eak interaction		| A or T/U  |
| M		| a**M**ino group		| A or C  |
| K		| **K**eto group		| G or T/U |   
| H		| not G		| A, C or T/U |   
| B		| not A		| C, G or T/U |  
| V		| not T/U		| A, C or G |  
| D		| not C		| A, G or T/U |  
| N		| a**N**y		| A, C, G or T/U |  
| - | none | Gap  

This table was adapted from [Johnson, 2010, *Bioinformatics*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/#B1).

An extended IUPAC genetic alphabet was also described in 2010 [(Johnson, 2010, *Bioinformatics*)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/#B1). The extended code uses additional characters, underlining, and bolding as well as the original 16 character code (all meanings maintained) to denote all possible combinations or relationships between bases. Among other uses, this has been valuable for representing genetic variation in DNA sequences. You can explore the details on the extended code in Tables 2 & 3 of [(Johnson, 2010)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865858/#B1).

If you're working with BioStrings objects, and need a reminder of the basic characters of the extended code, you can just type `IUPAC_CODE_MAP` when you have the BioStrings package loaded into your R session.

```r
# print iupac code to console
IUPAC_CODE_MAP
```

For now, we will use use the standard, unambiguous DNA bases (A, G, C, T) for some examples.
```{r}
# randomly sample from specific characters in DNA_ALPHABET to create a longer sequence
seq = sample(DNA_ALPHABET[c(1:4)], size=100, replace=TRUE)
seq

# use the paste command to collapse the characters into a single string
seq = paste(seq, collapse="")
seq

# use the DNAString constructor function to turn this sequence into a DNAString class object
seq.dnastring <- DNAString(seq)
seq.dnastring
```

Now collect some basic information on your sequence.
```r
# confirm how long it is
length(seq.dnastring)

# what is the frequency of each base in your sequence
alphabetFrequency(seq.dnastring, baseOnly=TRUE, as.prob=TRUE)

# what is the frequency of your favourite base (which is obviously Adenine)
letterFrequency(seq.dnastring, "A", as.prob=TRUE)

# return the frequency of dinucleotide pairs
dinucleotideFrequency(seq.dnastring, as.prob=TRUE)

# or even trinucleotides
trinucleotideFrequency(seq.dnastring, as.prob=TRUE)
```

We can also perform some basic manipulations of our sequence using BioStrings functions.
```r
# subset the sequence for 10 specific bases of interest  
seq.dnastring[10:19]

# this can also be done with the `subseq()` function  
subseq(seq.dnastring, 10, 19)

# get the reverse of our sequence
reverse(seq.dnastring)

# get the reverse COMPLEMENT sequence
reverseComplement(seq.dnastring)

# get the reverse complement of the first 10 bases in your sequence
reverseComplement(subseq(seq.dnastring, 1, 10))

# translate our DNA sequence
translate(seq.dnastring)
```

Again, our example is a little impractical since we are usually working with a set of sequences, for example the chromosomes in a reference genome. This is where the `DNAStringSet` object class becomes useful. `DNAStringSet` allows you to store, name, and manipulate multiple sequences in one BioStrings object.
```r
# remove the old single sequence from our global R environment
rm(seq)

# create a new variable, and fill it with individual sequences created as we did above
seq.dnass <- NULL
seq.dnass[1] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[2] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[3] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[4] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass[5] = paste(sample(DNA_ALPHABET[c(1:4)], size=50, replace=TRUE), collapse="")
seq.dnass

# how long is this object
length(seq.dnass)

# use the constructor function DNAStringSet to make the DNAStringSet object with your sequences
dna.st.set = DNAStringSet(seq.dnass)

# how long is the DNAStringSet object
length(seq.dnass)

# name all your sequences
names(seq.dnass) = paste("barcode-", 1:5, sep="")
seq.dnass
```

Like XString, a virtual class exists for `DNAStringSet` class objects called *XStringSet*, which also contains object classes for storing RNA and AA sequences (`RNAStringSet` and `AAStringSet`).

---

#### Working with *BSGenome* reference genomes

Now that we understand the major methods and classes implemented by BioStrings, lets load a complete reference genome and start exploring it.
```r
# assign the genome to a variable using getBSgenome() (you need to have the package for the BSgenome you are trying to load already installed)
genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
genome

# check the structure
str(genome)

# print the metadata for the genome
metadata(genome)
```

By default, the *BSGenomes* come with no sequence masking. It is common when woprking with reference genomes to mask regions that may contain ambiguous sequences, such as repeat regions, that you wish to ignore in your analyses. To obtain a masked genome, you should set `masked=TRUE` in the `getBSgenome()` function. This will load a genome in which specific sequences have been masked in a hierarchical fashion using the following criteria:  
1. Gaps in the genome assembly
2. Sequences with intra-contig ambiguities
3. regions flagged by [*RepeatMasker*](http://www.repeatmasker.org/)
4. regions flagged by [*Tandem Repeat Finder*](https://tandem.bu.edu/trf/trf.html)

Load the masked reference and compare to the unmasked version.
```r
genome.m <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")
class(genome.m)
genome.m

# unmaksed genome
class(genome)
genome
```

The masked genomes utilized the `MaskedXString` class implemented by *BioStrings* to denote the masked sequences. If you print a specific sequence to the console from the masked genome, you will also get a high-level summary of masking for that sequence.
```r
# return basic sequence information summary
seqinfo(genome.m)

# print chromosome 1
genome.m$chr1

# unmasked genome
seqinfo(genome)
genome$chr1
```

Let's move forward with the masked genome for today. Remove the `genome` variable from your working environment and replace it with the masked genome for convenience.
```r
rm(genome)

genome <- getBSgenome("BSgenome.Mmusculus.UCSC.mm10.masked")
```

We can perform all the same *BioStrings* based-methods on the sequences stored in our *BSGenome* object. For example:
```r
# assign chr 1
chr1 <- genome$chr1

# confirm how long it is
length(chr1)

# subset it
subseq(chr1, 1, 100)
subseq(chr1, 100498, 100598)

# what is the frequency of each base in your sequence
alphabetFrequency(chr1, baseOnly=TRUE, as.prob=TRUE)

# what is the frequency of your favourite base
letterFrequency(chr1, "A", as.prob=TRUE)
```

Beyond the basic BioStrings based methods, there is one important method implemented by the BSGenome using the `getSeq()` function, that extracts sequences at request from a BSgenome or XStringSet class object. We will use getSeq() functionality in our example below to demonstrate how you might use BSGenome packages in a typical NGS analysis.

---

#### Example: Extracting sequences flanking ChIP-seq peaks

Once peak regions have been identified to describe the potential binding locations of transcription factor (TF) or histone modification, a common task in the analysis of ChIP-seq data is to scan the sequences immediately surrounding these peaks in order to identify sequences enriched over these peak regions that may represent the binding motif for that TF. To achieve this, we need to obtain the sequences for these peaks from the reference genome that the samples were aligned to (mm10). The cartoon below depicts this overall workflow.

<p align="center">
<img src="../figures/motif-example.png" title="xxxx" alt="context"
	width="90%" height="90%" />
</p>

As an example, we will again use data from the ENCDOE project, where mouse forebrain tissues were ChIP'd for CTCF, a critical TF for diverse cellular processes that performs a wide range of transcriptional activation/repression functions at a genome-wide level. Called CTCF peaks for this experiment were downloaded from the ENCODE website [here](https://www.encodeproject.org/experiments/ENCSR677HXC/).

Read in the BED file as a *GRanges* object using *rtracklayer* function `import()` as we have done previously. We can then use the `getSeq()` function to return sequences from our previously assigned BSGenome object (UCSC - mm10, assigned to the variable *genome*) that cover the regions specified in the GRanges object.
```r
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
```

Since the object returned by `getSeq()` is a DNAStringSet class object, we can use BioStrings based methods to perform operations on the sequences directly. For example, we might be interested in checking the nucleotide frequencies across all peaks.
```r
# calculate nucleotide freqs.
nt_freqs <- alphabetFrequency(ctcf_seqs, baseOnly=TRUE, as.prob=TRUE)

# calculate mean value for nucleotide freqs across all peaks
round(apply(nt_freqs, 2, mean), digits=2)
```

We might also be interested in visualizing the distribution of the peak width, to get an idea of how much they vary. We can use the `width` accessor function to extract the width of each peak, and base R functions for plotting.
```r
hist(width(ctcf_seqs),
     col = "darkgray",
     xlab = "Peak width (bp)",
     main = "CTCF peak width distribution")
```

We could now export these sequences to a FASTA file (using `writeXStringSet()`) however several motif discovery softwares require that peaks be of the same size (width). To do this in a meaningful way for our ChIP-seq data, we will need to find the center of each peak, and then restrict to the a certain number of bases flanking either side of the center position. We will need to go back to the ranges from our original BED file to resize the peaks to the desired width around the center, then re-extract the sequences for those regions.
```r
# resize the regions from the BED file
bed_centered <- resize(bed, width = 400, fix = "center")
bed_centered

# check their with
width(bed_centered)

# extract sequences again
ctcf_seqs_cent <- getSeq(genome, bed_centered)
ctcf_seqs_cent
```

Now we are ready to export these sequences in FASTA file format, which is used as the default format as input to many motif discovery algorithms. As mentioned above, we can do this for DNAStringSet objects with the function `writeXStringSet()`.
```r
# add names to peaks in ctcf_seqs so that FASTA entries have names
names(ctcf_seqs) <- paste0(seqnames(bed), ":", start(bed), "-", end(bed))

# export peaks to FASTA file
writeXStringSet(ctcf_seqs, file="CTCF-peaks-resized.fa")
```

After you write the file, go to your the UNIX command line and have a look at your FASTA file to confirm it looks correct.

**Note:** There are other ways you could have performed this task outside of the functionality implemented by BioStrings and BSGenome. The major advantage of performing this analysis in R/Bioconductor is that you can leverage the methods and object classes implemented in BioStrings, BSGenome, GRanges, and other Bioconductor packages for your analysis.  

If you did not require R/Bioconductor functionality for further analysis, you could perform a similar analysis on the UNIX command line with [*bedtools*](https://bedtools.readthedocs.io/en/latest/) and its `getfasta` tool, which allows you to extract sequences from a BED/GTF/VCF file and export them to a FASTA file. Ofcourse, this do mean you would need to have a copy of the reference genome available to you.

---

#### Other functionality in BioStrings

BioStrings also provides functionality for a number of other analytical tasks that you may want to perform on a set of sequences stored using the *XString* and *XStringSet* method, for example:  
* trimming sequence ends based on pattern matching using `trimLRPatterns()`
* local and global alignment problems using `pairwiseAlignment()`
* read in multiple sequence alignments using `readDNAMultipleAlignment()`
* motif searches with a Position Weight Matrix (PWM) using `matchPWM()` (commonly done in ChIP-seq & ATAC-seq)
* palindrome searching using findPalindromes `findPalindromes()`
* computing edit distances between sets of sequences using `stringDist()`  

> ATTENTION: (performed at home or during extra time) that provides some exercises exploring additional functionality available in BioStrings can be found [here](https://github.com/Dartmouth-Data-Analytics-Core/Bioinformatics_workshop-Dec-2021/blob/master/Day-3/03-optional-exercise---BioStrings-extra-stuff.md).

An excellent BioStrings tutorial is available [here](https://bioconductor.org/help/course-materials/2011/BioC2011/LabStuff/BiostringsBSgenomeOverview.pdf) from one of the BioStrings creators, that covers much of the same material as we have above, but in more detail and with more complex examples.

> Similar tasks (and many other things) can be performed in python using tools like [*biopython*](https://biopython.org/). The major advantage of the *BioStrings* package is its interoperability with other R/BioConductor packages, such as the *BSGenome* package.
