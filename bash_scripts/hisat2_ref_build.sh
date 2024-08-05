#!/usr/bin/env bash 

## Root directory
# P_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/Pf3D7_genomes_gtfs_indexs
P_dir=/lustre/scratch126/tol/teams/lawniczak/users/jr35/genomes_gtfs_indexs

## Make relevant directories
mkdir -p ${P_dir}/Pf66_hisat_refs/genome_ref
mkdir -p ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref

## Copy fasta and GTF file used by sk21 to relevant location
mca_dir=/lustre/scratch126/tol/teams/lawniczak/projects/malaria_single_cell/mali_field_runs/2022/data/references_for_souporcell/
cp ${mca_dir}/PlasmoDB-66_Pfalciparum3D7_Genome.fasta ${P_dir}/Pf66_hisat_refs/genome_ref/genome.fasta
cp ${mca_dir}/PlasmoDB-66_Pfalciparum3D7_Genome.fasta ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.fasta

cp /lustre/scratch126/tol/teams/lawniczak/users/sk21/gafiles/db66/PlasmoDB-66_Pfalciparum3D7_proteincoding_ncRNAs.gtf ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.gtf

## Prepare genome fasta
hsat_dir=/software/team222/jr35/hisat2/hisat2/

${hsat_dir}/hisat2-build -p 16 ${P_dir}/Pf66_hisat_refs/genome_ref/genome.fasta ${P_dir}/Pf66_hisat_refs/genome_ref/genome

## Prepare transcript fasta
### Splice sites
${hsat_dir}/hisat2_extract_splice_sites.py ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.gtf > ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.ss

### Exons
${hsat_dir}/hisat2_extract_exons.py ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.gtf > ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.exon

### Build
${hsat_dir}/hisat2-build -p 16 --exon ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.exon --ss ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.ss ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran.fasta ${P_dir}/Pf66_hisat_refs/genome_w_tran_ref/genome_tran


##################### Pm

## Make relevant directories
mkdir -p ${P_dir}/Pm66_hisat_refs/genome_ref
mkdir -p ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref

## Copy fasta and GTF file used by sk21 to relevant location
cp ${mca_dir}/PlasmoDB-66_PmalariaeUG01_Genome.fasta ${P_dir}/Pm66_hisat_refs/genome_ref/genome.fasta
cp ${mca_dir}/PlasmoDB-66_PmalariaeUG01_Genome.fasta ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.fasta

cp /lustre/scratch126/tol/teams/lawniczak/users/sk21/gafiles/db66/PlasmoDB-66_PmalariaeUG01_proteincoding_ncRNAs.gtf ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.gtf

## Prepare genome fasta
hsat_dir=/software/team222/jr35/hisat2/hisat2/

${hsat_dir}/hisat2-build -p 16 ${P_dir}/Pm66_hisat_refs/genome_ref/genome.fasta ${P_dir}/Pm66_hisat_refs/genome_ref/genome

## Prepare transcript fasta
### Splice sites
${hsat_dir}/hisat2_extract_splice_sites.py ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.gtf > ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.ss

### Exons
${hsat_dir}/hisat2_extract_exons.py ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.gtf > ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.exon

### Build
${hsat_dir}/hisat2-build -p 16 --exon ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.exon --ss ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.ss ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran.fasta ${P_dir}/Pm66_hisat_refs/genome_w_tran_ref/genome_tran

