# Influenza A pipeline

Influenza A EQA Sequencing Pipeline
This is a bioinformatics pipeline for the analysis of Influenza A virus sequencing data, specifically designed for External Quality Assessment (EQA) purposes. The pipeline takes raw sequencing reads (Illumina or ONT) for H1N1 and H3N2 strains, performs quality control, alignment, variant calling, and generates a consensus sequence. Finally, it performs clade analysis and generates a comprehensive, strain-specific summary report.
The entire workflow is built using Nextflow, which ensures portability, reproducibility, and scalability.

Pipeline Overview
The pipeline performs the following steps:
Quality Control: Raw reads are assessed using FastQC.
Read Trimming/Filtering:
Illumina reads are trimmed for adapters and low-quality bases using Trimmomatic.
ONT reads are filtered for quality using NanoFilt.
Alignment: Reads are aligned to a strain-specific reference genome.
Illumina reads are aligned with BWA-MEM.
ONT reads are aligned with Minimap2.
Variant Calling: SNPs and Indels are called from the alignments using BCFtools.
Consensus Generation: A consensus FASTA sequence is generated for each sample based on the called variants.
Clade Analysis: The hemagglutinin (HA) segment is extracted from the consensus and analyzed with Nextclade to determine the genetic clade and identify key mutations.
Reporting: A separate, comprehensive summary report is generated for each strain (H1N1, H3N2) using MultiQC, aggregating the results from all previous steps.


Dependencies
The pipeline relies on Nextflow to run and Conda (specifically the mamba implementation for speed) to manage all software dependencies.

The main bioinformatics tools used are:
FastQC
Trimmomatic
NanoFilt
BWA
Minimap2
Samtools
BCFtools
Nextclade
MultiQC

Installation and Setup
1. Clone the Repository
git clone (https://github.com/ShabanaBi2025/neqas_pipeline)
cd [neqas_pipeline]

2. Create the Conda Environment
All the software required for the pipeline is defined in the envs/environment.yml file. Create the Conda environment using mamba:
mamba env create -f envs/environment.yml

3. Set Up the Directory Structure
The pipeline expects a specific directory structure for input files.
.
├── data/
│   ├── h1n1/
│   │   ├── sample1_R1.fastq.gz
│   │   └── sample1_R2.fastq.gz
│   └── h3n2/
│       └── sample2.fastq.gz
├── references/
│   ├── h1n1/
│   │   └── h1n1_reference.fasta
│   └── h3n2/
│       └── h3n2_reference.fasta
├── adapters/
│   └── TruSeq3-PE.fa           # Adapter file for Trimmomatic
├── envs/
│   ├── environment.yml
├── modules/
│   ├── qc.nf
│   ├── trimmomatic.nf         
│   ├── nanofilt.nf             
│   ├── align_illumina.nf
│   ├── align_ont.nf
│   ├── mutltiqc.nf         
│   ├── bcftools.nf             
│   └── nextclade.nf
├── main.nf
└── nextflow.config             # Main config includes 

5. Prepare Reference Genomes (One-time step)
Before running the main pipeline, you must clean and index your reference genomes. This is a one-time setup step.

nextflow run main_clean_refs.nf

This will create a results/cleaned_reference directory containing the prepared genomes that the main pipeline will use.

5. Download Nextclade Datasets (One-time step)
You also need to download the Nextclade datasets for H1N1 and H3N2.

# Download H3N2 dataset
nextclade dataset get --name 'flu_h3n2_ha' --output-dir 'nextclade_datasets/h3n2'

# Download H1N1 dataset
nextclade dataset get --name 'flu_h1n1pdm_ha' --output-dir 'nextclade_datasets/h1n1'

Configuration
All pipeline parameters are defined in the nextflow.config file. The most important parameters are set by default but can be modified if needed:

params.raw_reads: The path pattern to find the input FASTQ files (default: data/**/*.fastq*).

params.outdir: The directory where all results will be saved (default: results).

Usage
Once the setup is complete, you can run the entire pipeline with a single command from the project's root directory.

nextflow run main.nf

To re-run the pipeline from the last failed step, use the -resume flag:

nextflow run main.nf -resume

Output
The pipeline will create a results/ directory containing the outputs from each step in its own sub-folder (e.g., aligned/, variants/, nextclade/).

The final, aggregated reports can be found in:
results/multiqc_report/

This directory will contain a separate folder for each strain (e.g., h1n1/, h3n2/), and inside each, you will find a multiqc_report.html file. This HTML file can be opened in any web browser and provides a comprehensive, interactive summary of all the results for that specific strain.
