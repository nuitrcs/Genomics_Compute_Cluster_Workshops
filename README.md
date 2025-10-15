# Genomics Compute Cluster Workshops

As part of the Genomics Compute Cluster (GCC) services, Northwestern IT's Research Computing and Data Services (RCDS) offers workshops tailored for the GCC community.
Materials for these workshops are available to the public and can be found linked from this page.

## Foundations in Genomics Analyses

This workshops series has been taught in the Fall quarters of 2024, and 2025. Topics vary slightly between years but generally include basic skills for performing computational analysis of sequence data on Northwestern's HPC system, Quest.

### Getting Genomic Data onto Quest
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
  Public genomic repositories such as the Gene Expression Omnibus (GEO), Sequence Read Archive (SRA), and the European Nucleotide Archive (ENA) are invaluable resources but can be challenging to use due to their diverse structures, metadata formats, and download protocols. This workshop introduces key tools and workflows for accessing raw sequencing data from major public databases to work with on Quest. Participants will learn how to map between GEO and SRA accessions and retrieve metadata and sequence files using command-line tools. We will also discuss best practices for data management, metadata parsing, and troubleshooting common issues in retrieval workflows.
  [Materials available here.](https://github.com/palupaca/Getting-Genomic-Data-Onto-Quest)
</details>

### Sequencing Filetypes and Quality Control
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
  Sequencing data is essential to many computational biology studies. However, how to use these data types is often not covered in standard curriculum. This workshop will cover the structure and use of common sequencing filetypes, including the difference between fasta and fastq, and zipped and unzipped files. We will also work through calculating quality control statistics on these filetypes with fastqc.
  [Materials available here.](https://github.com/nuitrcs/genomic_filetypes)
</details>

### Sequence Alignment and Mapping DNA-seq Reads
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
The first step in almost all bioinformatic pipelines is sequence alignment. There are many software tools available to accomplish this task. Choosing an appropriate tool and then using it at scale can be intimidating. This workshop will remove some of the mystery by suggesting aligners for specific tasks and covering the basics of using each suggested tool including: bwa, bowtie2, and minimap2. Participants will compose and run submission scripts for each of these, and we will look at the output they create.
  [Materials available here.](https://github.com/nuitrcs/sequence_alignment)
</details>
### RNA-seq Alignment
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
Read alignment or mapping is a computational process to determine where in the reference genome the short RNA reads originated from. This workshop will introduce how to perform read alignment using STAR (Spliced Transcripts Alignment to a Reference), a splice-aware aligner designed to specifically tackle many challenges involved in RNA-seq data mapping. We will introduce STAR’s alignment strategy and interactively work through each step of alignment using an example dataset. We will also discuss how STAR is incorporated in many Nextflow pipelines (for example nf-core/rnaseq).
  [Materials available here.]()
</details>

### Genome Browsers
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
This workshop provides a brief overview of the most popular features of two genome browsers: the UCSC Genome Browser and the Ensembl Genome Browser. Both browsers allow the user to retrieve comprehensive information on a gene or genomic sequence in context of the genome, with multiple customizable tracks to display curated data from numerous external sources for features such as polymorphisms, transcript variants and histone modifications.
This workshop will cover how to locate these browsers, manage tracks and annotations, get sequence data from the browser view, and download data from each site. We will compare how each browser displays data, so the user can decide which is best for their purposes.
 [Materials available here.]()
</details>

### Sequence Similarity Searching
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
Sequence similarity searches can be done in multiple ways on multiple platforms. It is useful for comparing or discovering conserved regions across sometimes very dissimilar sequences. BLAST (Basic Local Alignment Search Tool) aligns a sequence (nucleotide or peptide) to a database of other sequences, or can align two sequences to each other. This workshop will cover the basics of sequence similarity searching with NCBI’s BLAST, as well as introduce specialized BLAST tools that can help you find statistically significant matches to a nucleotide or protein query sequence, discover homology across species, target your search to specific taxonomic groups in the BLAST database, and retrieve data for further analysis.
 [Materials available here.]()
</details>

### Setting up an R Environment for Analysis with the GCC
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
R is a popular language for data analysis and many R packages have been written for genomic analysis, especially analysis of single cell RNA-seq data. Installing and managing these packages is relatively straightforward on a local laptop but performing these analyses with large genomic or transcriptomic datasets often requires more computational power than an individual laptop provides. This workshop will cover R package installation on Quest and interfacing with a custom R environment via RStudio on Quest OnDemand. We will focus on installation of Bioconductor packages and packages used with scRNA-seq data.
 [Materials available here.]()
</details>

### Getting Started with nf-core Nextflow Pipelines
<details markdown="1">
  <summary>Abstract and Link to Materials</summary>
Nextflow is a workflow management tool designed to enable the creation of reproducible workflows that function across computational resources through the use of software containers. Many curated bioinformatics Nextflow pipelines can be found through the nf-core. In this workshop, we will explore their collection and work through running an nf-core pipeline on Quest. The biggest consideration for using Nextflow on Quest is setting the configuration to use Quest’s resources appropriately. We will cover how this is done through the nu-genomics profile for nf-core pipelines as well as writing a custom configuration file when necessary.
 [Materials available here.]()
</details>

## Topics in Computational Genomics
