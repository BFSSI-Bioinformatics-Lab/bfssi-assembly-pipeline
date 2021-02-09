#!/usr/bin/env nextflow

// Default parameters
params.outdir = './pipeline_results'
params.reads = './reads/*_{R1,R2}.fastq.gz'

// QC params
params.fastqc = false
params.quast = false

// Flags to run optional post-processing steps
params.plasmids = false
params.staramr = false
params.annotate = false
params.mlst = false
params.abricate = false
params.abricate_db = 'ncbi'

reads = Channel
    .fromFilePairs( params.reads , flat:true)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch }

reads = Channel
    .fromFilePairs( params.reads , flat:true)
    .set { reads_prefilter_ch }

process FastQCPre {
    container 'staphb/fastqc:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'

    when:
    params.fastqc

    input: 
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from reads_prefilter_ch

    output:
    path('fastqc/prefilter/*')

    script:
    """
    mkdir -p ./fastqc/prefilter
    fastqc -o ./fastqc/prefilter $fwd_reads $rev_reads 
    """
}

process runBBMapReadRepair {
    container 'staphb/bbtools:latest'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from read_pairs_ch

    output:
    tuple val(pair_id), path('*_R1.repaired.fastq.gz'), path('*_R2.repaired.fastq.gz') into read_pairs_repaired_ch

    script:
    """
    repair.sh in=$fwd_reads in2=$rev_reads out1=${pair_id}_R1.repaired.fastq.gz out2=${pair_id}_R2.repaired.fastq.gz overwrite=t
    """
}


process runBBMapAdapterTrimming {
    container 'staphb/bbtools:latest'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from read_pairs_repaired_ch

    output:
    tuple val(pair_id), path('*_R1.trimmed.fastq.gz'), path('*_R2.trimmed.fastq.gz') into read_pairs_trimmed_ch

    script:
    """
    bbduk.sh in1=$fwd_reads in2=$rev_reads out1=${pair_id}_R1.trimmed.fastq.gz out2=${pair_id}_R2.trimmed.fastq.gz  ref=adapters tpe tbo overwrite=-t unbgzip=f ktrim=r k=23 mink=1 hdist=1
    """
}

process runBBMapQualityFiltering {
    container 'staphb/bbtools:latest'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from read_pairs_trimmed_ch

    output:
    tuple val(pair_id), path('*_R1.qfilter.fastq.gz'), path('*_R2.qfilter.fastq.gz') into read_pairs_qfilter_ch

    script:
    """
    bbduk.sh in1=$fwd_reads in2=$rev_reads out1=${pair_id}_R1.qfilter.fastq.gz out2=${pair_id}_R2.qfilter.fastq.gz unbgzip=f qtrim=rl trimq=10
    """
}

process runBBMapReadCorrection {
    container 'staphb/bbtools:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'

    input:
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from read_pairs_qfilter_ch

    output:
    tuple val(pair_id), path('*_R1.corrected.fastq.gz'), path('*_R2.corrected.fastq.gz') into read_pairs_corrected_ch, reads_postqc_ch

    script:
    """
    tadpole.sh in1=$fwd_reads in2=$rev_reads out1=${pair_id}_R1.corrected.fastq.gz out2=${pair_id}_R2.corrected.fastq.gz mode=correct
    """
}

process FastQCPost {
    container 'staphb/fastqc:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'

    when:
    params.fastqc

    input: 
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from reads_postqc_ch

    output:
    path('fastqc/postfilter/*')

    script:
    """
    mkdir -p ./fastqc/postfilter
    fastqc -o ./fastqc/postfilter $fwd_reads $rev_reads 
    """
}

process runSKESA {
    container 'staphb/skesa:latest'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fwd_reads), path(rev_reads) from read_pairs_corrected_ch

    output:
    tuple val(pair_id), path(fwd_reads), path(rev_reads), path('*.fasta') into assembly_ch

    script:
    """
    skesa --use_paired_ends --fastq "$fwd_reads,$rev_reads" --contigs_out ${pair_id}.fasta
    """
}

process produceBAM {
    container 'staphb/bbtools:latest'
    tag "$pair_id"

    input:
    tuple val(pair_id), path(fwd_reads), path(rev_reads), path(assembly) from assembly_ch

    output:
    tuple val(pair_id), path(assembly), path('*.bam') into bam_ch

    script:
    """
    bbmap.sh in1=$fwd_reads in2=$rev_reads ref=$assembly out=${pair_id}.bam overwrite=t deterministic=t
    """
}

process sortBAM {
    container 'staphb/samtools:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'

    input:
    tuple val(pair_id), path(assembly), path(bamfile) from bam_ch

    output:
    tuple val(pair_id), path(assembly), path('*.sorted.bam') into sorted_bam_ch

    script:
    """
    samtools sort --output-fmt BAM -o ${pair_id}.sorted.bam $bamfile
    """
}

process indexBAM {
    container 'staphb/samtools:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'


    input:
    tuple val(pair_id), path(assembly), path(sorted_bamfile) from sorted_bam_ch

    output:
    tuple val(pair_id), path(assembly), path(sorted_bamfile), path('*.bai') into indexed_bam_ch, bamqc_ch

    script:
    """
    samtools index $sorted_bamfile
    """
}

process runPilon {
    container 'staphb/pilon:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'

    input:
    tuple val(pair_id), path(assembly), path(indexed_bamfile), path(bam_indexfile) from indexed_bam_ch

    // Output final polished assembly into n separate channels for post-processing steps
    output:
    tuple val(pair_id), path(assembly) into polished_assembly_ch1,
                polished_assembly_ch2,
                polished_assembly_ch3,
                polished_assembly_ch4,
                polished_assembly_ch5,
                polished_assembly_ch6,
                polished_assembly_ch7

    script:
    """
    pilon --genome $assembly --bam $indexed_bamfile --outdir $params.outdir --output $pair_id
    """
}

process runQuast {
    container 'staphb/quast:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id", mode: 'symlink'

    when:
    params.quast

    input:
    tuple val(pair_id), path(assembly) from polished_assembly_ch7

    output:
    path('quast/*')

    script:
    """
    mkdir -p ./quast
    quast.py -o ./quast --labels $pair_id --glimmer --no-plots --no-icarus --no-snps $assembly
    """

}

process runProkka {
    container 'staphb/prokka:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id/", mode: 'symlink'

    when:
    params.annotate

    input:
    tuple val(pair_id), path(assembly) from polished_assembly_ch1

    output:
    path('prokka/*')

    script:
    """
    mkdir -p prokka
    prokka --centre BFSSI --compliant --prefix $pair_id --locustag $pair_id --force --outdir ./prokka $assembly
    """
}

process runStarAMR {
    container 'staphb/staramr:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id/", mode: 'symlink'

    when:
    params.staramr

    input:
    tuple pair_id, file(assembly) from polished_assembly_ch2

    output:
    path('staramr/*')

    // StarAMR currently hard coded to be very lenient with possible input genomes, maybe expose some parameters in the future
    // StarAMR can't deal with an existing directory, so -resume will break. Have to attempt to remove dir first.
    script:
    """
    rm -rf ./staramr
	staramr search --genome-size-lower-bound 100 --genome-size-upper-bound 8000000 --minimum-N50-value 5000 --minimum-contig-length 250 -o ./staramr $assembly
    """
}


process runAbricate {
    container 'staphb/abricate:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id/", mode: 'symlink'

    when:
    params.abricate

    input:
    tuple pair_id, file(assembly) from polished_assembly_ch6
    val(abricate_db) from params.abricate_db

    output:
    path('abricate/*')

    // Abricate supports numerous databases; allow user to specify which one they'd like (default ncbi)
    script:
    """
    mkdir -p abricate
    abricate --db $abricate_db $assembly > ./abricate/${pair_id}_${abricate_db}.tab
    """
}

process runMobRecon {
    container 'kbessonov/mob_suite:3.0.1'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id/", mode: 'symlink'

    when:
    params.plasmids

    input:
    tuple pair_id, file(assembly) from polished_assembly_ch3

    output:
    path('plasmids/*')

    script:
    """
    mkdir -p plasmids
    mob_recon --infile $assembly -o ./plasmids --run_typer --force
    """
}

process runMLST {
    container 'staphb/mlst:latest'
    tag "$pair_id"
    publishDir "$params.outdir/$pair_id/mlst", mode: 'symlink'

    when:
    params.mlst

    input:
    tuple pair_id, file(assembly) from polished_assembly_ch4

    output:
    path('*_mlst.tsv')

    script:
    """
    mlst $assembly > ${pair_id}_mlst.tsv
    """
}
