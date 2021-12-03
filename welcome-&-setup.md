# Welcome to the DAC Fundamentals of Bioinformatics workshop #

Before you attend the workshop there are a couple of things we would like you to do to get setup for using the tools that will be required during the course of the workshop. Please read through each of the sections below to ensure you are prepared to attend the workshop. We strongly recommend that you run through these steps several days in advance of the workshop, in case any troubleshooting is required.

---

## The terminal emulator ##

A terminal emulator is a more streamlined way of navigating a computational environment (once you get the hang of it). We will cover some basic commands to help orient you with using the terminal to interact with your machine on Day 1 of the workshop.

If you are using a Mac there is a terminal emulator program pre-installed on your machine, if you are using another OS we have included some links to popular terminal emulator programs below. Please select one of them download the program and open it up for the next step.

Operating system| Terminal emulators
---|---
Mac| Terminal (comes pre-installed)
Windows| [MobaXterm](https://mobaxterm.mobatek.net/download.html) <br> [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)
Linux| Konsole, Terminal, etc. (should be pre-installed but depends on the desktop environment you are running)

---

## The discovery HPC system ##

For those of you that indicated that you did not have an account on *discovery* you will need to request one [here](https://rc.dartmouth.edu) **AT LEAST** two (business) days before the workshop begins. There is a green button at the top of the page that says **Request an Account**, once you click the button you will be prompted to log in with your netID and password and then you can fill in the form to request your account.

Once you have a discovery account you can follow along with this video [here](https://youtu.be/VoHBlblsQfg) to log onto discovery using the command line.

To log onto discovery we will use the secure shell command `ssh`. 

```bash
ssh netID@discovery.dartmouth.edu
```

You will be prompted for a password and when you start typing nothing will show up in the prompt, I assure you though your keystrokes are recorded and you will be logged onto the discovery HPC environment if your password is correct. If your password is incorrect you will receive the warning "Permission denied, please try again." and will be prompted to enter your password again.

---

## Setting up a Conda Environment ## 

Once you are logged onto discovery you can load all of the software we will need for the workshop. 

Conda is a package management system that helps you find, install, and organize groups of packages needed for a particular task. Conda environments are really useful when working on HPC environments like Dartmouth's Discovery system because you can install packages locally without needing administrator permission. Conda environments are also useful for project continuity, the versions of the packages that you install in your environment and all of their dependencies will remain the same (unless you update them). We will be using a conda environment to make sure we all have the same version of many different bioinformatics software programs available to us. 

Before you begin using conda environments on discovery you will need to ensure that you have run the source code to enable the conda commands. Ensure you are logged into discovery and run the following command:

```bash
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
```

We recommend that you add the above line of code to your `.bashrc` file in your home directory, otherwise you will need to run this command each time you start a new session on discovery. To do this use the `nano` text editing program on discovery (this comes pre-installed) to copy the line above into your `.bashrc` file (we will talk more about what this file is and how to use it on Day 1).

```bash
# open the file with the Nano text editor
nano .bashrc

# copy this line to the file : source /optnfs/common/miniconda3/etc/profile.d/conda.sh

# use ctrl+x to exit the editor, then type "Y" to save the changes you made, then press Enter write the changes to .bashrc
```


Next you will have to run the following command to create a .conda/ directory in your home directory. This directory will store all of your personal conda environments, including the one we are about to build for this workshop. **You only have to run this command once to make this directory, so it does not need to be added to your .bashrc file.**

```bash
cd ~
mkdir -p .conda/pkgs/cache .conda/envs
```

Lastly you will need to create the conda environment that we will be using for the course in your personal directory of accessible conda environments. This takes about 15 minutes to execute and you will see all of the packages that are loaded into this environment. The number of packages being installed should indicate why conda environments are so useful, imagine having to load all of these packages individually it is much easier to load them with a single command in a conda environment.

```bash
conda env create -f /dartfs-hpc/scratch/fund_of_bioinfo/environment.yml
```

When you are ready activate the conda environment, which you will need for the work we are doing for days 1 and 2 of the workshop, you can use the following command. 

```bash
conda activate bioinfo
```

You will see that the activate command has worked when it reads (bioinfo) rather than (base) to the left of the prompt. When you are finished using a conda environment it is good practice to deactivate your session with the following command.

```bash
conda deactivate
```

You can follow along with me as I run these commands to create a new conda environment in this [video](https://youtu.be/73IHZlFvb5Q)

---

## Installing an SFTP client ##

**This is optional** but for those of you that are new to the command line this might be an easier way to move files between the HPC environment and your local machine. An SFTP client stands for secure file transfer protocol and will enable you to drag and drop files as you might in a finder window between your local machine and a remote location. I use FileZilla, which I believe works on Mac, Windows, and linux operating systems. You can download [FileZilla](https://filezilla-project.org/download.php?show_all=1) by following the link and selecting the version that is correct for your OS, then open the program to ensure that you have downloaded it successfully. Once you have Filezilla installed you can use this [video](https://youtu.be/A8w8Uw1OILA) to guide you in linking the SFTP client to your account on discovery.

---

## Install the Integrative Genomics Viewer (IGV)

We will be using the [Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/), a genome browser produced by researchers at the Broad Institute, to explore and visualize genomics data. 

<img src="figures/igv.png" height="100" width="100"/>

You will need to download and install the IGV desktop application for your operating system before the workshop begins. The latest versions of IGV can be found at their [downloads page](http://software.broadinstitute.org/software/igv/download). After installing IGV, try opening the application on your computer to confirm the installation was successful. 

---

## Setting up an R project ##

We will be using R-Studio to explore and analyze genomics data on day 2 and 3, therefore we ask that you have R and R-Studio installed prior to attending the workshop. You will need to be running at least version 3.6 to ensure all of the packages needed will run smoothly. The latest versions for R and R-Studio can be found [here](https://cran.r-project.org) and [here](https://rstudio.com/products/rstudio/download/).

Next you will need to set up a new project in R-Studio for this workshop. Projects in R are like containers for various jobs that you will perform, the history of the project will be loaded when you open a new project. By using a project you can install all of the tools ahead of time and they will be there for you when you need them during the workshop. In R-Studio under File select New directory and then select New Project and name the project something you will remember (bioinfo_workshop).

Now that you have your project loaded, run the following code to install all of the packages we will be using during the workshop. For those of you that haven't used RStudio before we have made a [video](https://youtu.be/UtZHS-q7buI) showing the successful installation of the R packages you will need using the commands below. 

I bumbled the description of the code chunks with the nested loops, so here is a better description for those that are interested. There are two loops and each loop starts with an if statement. The first loop states "if `biomaRt` is not installed enter this loop" and the second one "if `BioCManager` is not installed enter this loop", when the condition is fulfilled (the package is *not* installed) the loop is entered and the function `install.packages` is used to install the package. Each loop is exited once the packages are installed and the package is loaded with the 'library' function to make the functions contained in the package available during the current session.


```r
if (!any(rownames(installed.packages()) == "biomaRt")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(biomaRt)

if (!any(rownames(installed.packages()) == "IRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("IRanges")
}
library(IRanges)

if (!any(rownames(installed.packages()) == "GenomicRanges")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicRanges")
}
library(GenomicRanges)

if (!any(rownames(installed.packages()) == "Gviz")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("Gviz")
}
library(Gviz)

if (!any(rownames(installed.packages()) == "org.Hs.eg.db")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)

if (!any(rownames(installed.packages()) == "EnsDb.Hsapiens.v86")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("EnsDb.Hsapiens.v86")
}
library(EnsDb.Hsapiens.v86)

if (!any(rownames(installed.packages()) == "GenomicFeatures")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("GenomicFeatures")
}
library(GenomicFeatures)

if (!any(rownames(installed.packages()) == "VariantAnnotation")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("VariantAnnotation")
}
library(VariantAnnotation)

if (!any(rownames(installed.packages()) == "TxDb.Hsapiens.UCSC.hg38.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
}
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

if (!any(rownames(installed.packages()) == "TxDb.Mmusculus.UCSC.mm10.knownGene")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "BSgenome")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("BSgenome")
}
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

if (!any(rownames(installed.packages()) == "ChIPseeker")){
  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("ChIPseeker")
}
library(ChIPseeker)

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

BiocManager::install("BSgenome.Mmusculus.UCSC.mm10.masked")

sessionInfo()
```

---

## Downloading the data ##

The commands that you will be following can be found in markdown `(.md)` files where there is a brief description of each command and how it is applied to the data and what it does followed by an example command that you can copy and paste into the terminal window. The majority of day 1 and 2 will be using the terminal window on your local machine, with an open `ssh` connection to discovery, as we will be running `bash` code. For some of day 2 and most of day 3 you will be using RStudio on your local machine to run the commands in the markdown files (`.md`) located in this GitHub repo. 


In your terminal window **on your local machine** navigate to where you want to download the files needed for this workshop. 

**On Monday** before you log onto the first zoom session we will make the workshop materials public and you should download those to your local machine (preferably in the same location as you downloaded the setup materials) with the following command: 

```bash
git clone https://github.com/Dartmouth-Data-Analytics-Core/Bioinformatics_workshop-Dec-2021/
```

---

If you have issues with any part of the installation and setup, please reach out to us directly (contact details are in the workshop ReadMe page) or come to bioinformatics help hours **December 3, 2021 12:30-1:30PM** the link is in your welcome email. 

