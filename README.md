# Genomics Compute Cluster Workshops

As part of the Genomics Compute Cluster (GCC) services, Northwestern IT's Research Computing and Data Services (RCDS) offers workshops tailored for the GCC community.
Materials for these workshops are available to the public and can be found linked from this page.

## Foundations in Genomics Analyses

This workshops series has been taught in the Fall quarters of 2024, and 2025. Topics vary slightly between years but generally include basic skills for performing computational analysis of sequence data on Northwestern's HPC system, Quest.

<details markdown="1" open>
  <summary markdown="1">
    
### Command Line Introduction
  </summary>
  Working from the command line gives computational biologists access to many analysis tools and enables the creation of reproducible scripts. This workshop will introduce commands for file system navigation, teach attendees how to create and edit scripts from the terminal with the text editor nano, and introduce utilities for interacting with our high-performance computer cluster, Quest.
  
  [Materials available here.](https://github.com/nuitrcs/gcc_command_line)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Getting Genomic Data onto Quest
  </summary>
  Public genomic repositories such as the Gene Expression Omnibus (GEO), Sequence Read Archive (SRA), and the European Nucleotide Archive (ENA) are invaluable resources but can be challenging to use due to their diverse structures, metadata formats, and download protocols. This workshop introduces key tools and workflows for accessing raw sequencing data from major public databases to work with on Quest. Participants will learn how to map between GEO and SRA accessions and retrieve metadata and sequence files using command-line tools. We will also discuss best practices for data management, metadata parsing, and troubleshooting common issues in retrieval workflows.
  
  [Materials available here.](https://github.com/palupaca/Getting-Genomic-Data-Onto-Quest)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Sequencing Filetypes and Quality Control
  </summary>
    Sequencing data is essential to many computational biology studies. However, how to use these data types is often not covered in standard curriculum. This workshop will cover the structure and use of common sequencing filetypes, including the difference between fasta and fastq, and zipped and unzipped files. We will also work through calculating quality control statistics on these filetypes with fastqc.
  
  [Materials available here.](https://github.com/nuitrcs/genomic_filetypes)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Sequence Alignment and Mapping DNA-seq Reads
    
  </summary>
The first step in almost all bioinformatic pipelines is sequence alignment. There are many software tools available to accomplish this task. Choosing an appropriate tool and then using it at scale can be intimidating. This workshop will remove some of the mystery by suggesting aligners for specific tasks and covering the basics of using each suggested tool including: bwa, bowtie2, and minimap2. Participants will compose and run submission scripts for each of these, and we will look at the output they create.
  
  [Materials available here.](https://github.com/nuitrcs/sequence_alignment)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### RNA-seq Alignment
    
  </summary>
    
Read alignment or mapping is a computational process to determine where in the reference genome the short RNA reads originated from. This workshop will introduce how to perform read alignment using STAR (Spliced Transcripts Alignment to a Reference), a splice-aware aligner designed to specifically tackle many challenges involved in RNA-seq data mapping. We will introduce STAR’s alignment strategy and interactively work through each step of alignment using an example dataset. We will also discuss how STAR is incorporated in many Nextflow pipelines (for example nf-core/rnaseq).
  
  [Materials available here.](https://github.com/nuitrcs/star_aligner_workshop)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Genome Browsers
    
  </summary>
    
This workshop provides a brief overview of the most popular features of two genome browsers: the UCSC Genome Browser and the Ensembl Genome Browser. Both browsers allow the user to retrieve comprehensive information on a gene or genomic sequence in context of the genome, with multiple customizable tracks to display curated data from numerous external sources for features such as polymorphisms, transcript variants and histone modifications.
This workshop will cover how to locate these browsers, manage tracks and annotations, get sequence data from the browser view, and download data from each site. We will compare how each browser displays data, so the user can decide which is best for their purposes.
 
  [Materials available here.](https://github.com/galterdatalab/foundations-genomebrowsers)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Sequence Similarity Searching
    
  </summary>
    
Sequence similarity searches can be done in multiple ways on multiple platforms. It is useful for comparing or discovering conserved regions across sometimes very dissimilar sequences. BLAST (Basic Local Alignment Search Tool) aligns a sequence (nucleotide or peptide) to a database of other sequences, or can align two sequences to each other. This workshop will cover the basics of sequence similarity searching with NCBI’s BLAST, as well as introduce specialized BLAST tools that can help you find statistically significant matches to a nucleotide or protein query sequence, discover homology across species, target your search to specific taxonomic groups in the BLAST database, and retrieve data for further analysis.
 
  [Materials available here.](https://github.com/galterdatalab/foundations-sequence-similarity)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Setting up an R Environment for Analysis with the GCC
    
  </summary>
    
R is a popular language for data analysis and many R packages have been written for genomic analysis, especially analysis of single cell RNA-seq data. Installing and managing these packages is relatively straightforward on a local laptop but performing these analyses with large genomic or transcriptomic datasets often requires more computational power than an individual laptop provides. This workshop will cover R package installation on Quest and interfacing with a custom R environment via RStudio on Quest OnDemand. We will focus on installation of Bioconductor packages and packages used with scRNA-seq data.
 
  [Materials available here.](https://github.com/nuitrcs/R_environments_GCC)
  
</details>


<details markdown="1" open>
  <summary markdown="1">
    
### Getting Started with nf-core Nextflow Pipelines
    
  </summary>
    
Nextflow is a workflow management tool designed to enable the creation of reproducible workflows that function across computational resources through the use of software containers. Many curated bioinformatics Nextflow pipelines can be found through the nf-core. In this workshop, we will explore their collection and work through running an nf-core pipeline on Quest. The biggest consideration for using Nextflow on Quest is setting the configuration to use Quest’s resources appropriately. We will cover how this is done through the nu-genomics profile for nf-core pipelines as well as writing a custom configuration file when necessary.
 
  [Materials available here.](https://github.com/nuitrcs/nextflow_nfcore_intro)
  
</details>

## Topics in Computational Genomics

<details markdown="1" open>
  <summary markdown="1">
    
### Peak Calling with MACS2
    
  </summary>
    
Peak calling is a computational method used to identify areas in the genome that have been enriched with aligned reads in ChIPseq and ATACseq datasets. This workshop will introduce Model-based Analysis of ChIP-Seq (MACS2), one of the commonly used peak callers, and demonstrate how to call peaks from aligned reads on Quest. We will introduce some basic MACS2 parameters and interactively work through each step of peak calling using an example dataset.  

  [Materials available here.](https://github.com/nuitrcs/MACS2_workshop)
  
</details>

<details markdown="1" open>
  <summary markdown="1">
    
### Scaling Up for High-Throughput Computing
    
  </summary>
    
Computational genomics workflows are most useful when they handle many samples. This workshop will introduce bash variables and job arrays to help you scale up as efficiently as possible with the SLURM scheduling software on Quest. We will interactively work through converting scripts that analyze one sample to analyze any number of samples, so you can launch one script that will handle all your sequencing files.

  [Materials available here.](https://github.com/nuitrcs/high_throughput_computing)
  
</details>

<details markdown="1" open>
  <summary markdown="1">
    
### Virtual Environments for Single-Cell Analysis
    
  </summary>
    
The number of software packages designed for single-cell RNA-seq analysis has increased dramatically since the advent of both the sequencing technology and more widely available GPUs. Setting up an environment with the correct software to use the GPUs on Quest requires attention to specific versions. In this workshop, we will walk through setting up a software environment with mamba that can be used on our GPUs from both the command line and from JupyterLab via Quest OnDemand, and we will discuss software version control more generally. 

  [Materials available here.](https://github.com/nuitrcs/virtual_environments_for_single_cell_analysis)
  
</details>

<details markdown="1" open>
  <summary markdown="1">
    
### Spatial Transcriptomics with Scanpy and Squidpy
    
  </summary>
    
Scanpy (single-cell analysis in Python) is a toolkit for single-cell gene expression data analysis using the programming language Python. This workshop will cover a basic tutorial for the use of Scanpy and related packages for analysis and visualization of spatial single-cell gene expression data, with a focus on explaining data types, and how to find and adjust different parameters of the included functions. 

  [Materials available here.](https://github.com/nuitrcs/spatial_transcriptomics)
  
</details>

