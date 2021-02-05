# BFSSI Assembly Pipeline

This is a prokaryotic assembly pipeline used internally at the BFSSI lab at Health Canada.
It is intended to be a simple, configurable Nextflow pipeline that will produce high
quality prokaryotic assemblies with useful post-processing.

The pipeline currently only expects paired-end reads, though hybrid assembly support is on the TODO list.

Currently, the Docker images used in this pipeline are sourced from
[StaPH-B](https://hub.docker.com/r/staphb/). By using Docker images for each container step, 
this pipeline allows the user to avoid onerous bioinformatics software 
installation and configuration. 

## Pipeline Overview

1. QC on reads with [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) (adapter trimming/quality filtering)
2. Error-correction of reads with [Tadpole](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/tadpole-guide/)
3. Assembly of short-reads with [SKESA](https://github.com/ncbi/SKESA)
4. Alignment of error-corrected reads against draft assembly with [BBmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
5. Polishing of assembly with [Pilon](https://github.com/broadinstitute/pilon/wiki)

### Optional post-processing
- [MLST profiling](https://github.com/tseemann/mlst) `--mlst`
- [Gene annotation](https://github.com/tseemann/prokka) `--annotate`
- [Plasmid detection](https://github.com/phac-nml/mob-suite) `--plasmids`
- [AMR profiling](https://github.com/phac-nml/staramr) `--amr`

## Installation/Requirements
Running this pipeline requires only `Nextflow` and `Docker`!

- [Nextflow installation instructions](https://www.nextflow.io/)
- [Docker installation instructions](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)

## Usage

Running the pipeline on paired-end reads with MLST and Prokka options enabled
`nextflow pipeline.nf --outdir ./results --reads "./data/*_{R1,R2}.fastq.gz" --mlst --annotate `

## TODO:
- Add pre- and post-assembly QC (FastQC, Quast, Qualimap)
- Add full support for plasmid/AMR detection
- Add support for long-read/hybrid assemblies (first with Unicycler/Trycycler)
- Provide a nice, interactive summary report