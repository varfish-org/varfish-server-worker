#!/usr/bin/env bash

set -euo pipefail

export DATA_DIR=${DATA_DIR-/tmp/mehari-data}
export SEQREPO_ROOT_DIR=$DATA_DIR/seqrepo

set -x

mkdir -p $DATA_DIR/tmp/$GENOME_RELEASE
cd $DATA_DIR/tmp/$GENOME_RELEASE

wget -q \
     https://ftp.ensembl.org/pub/$ENSEMBL_RELEASE/fasta/homo_sapiens/cdna/Homo_sapiens.$ENSEMBL_TOKEN.cdna.all.fa.gz
wget -q \
    https://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.{1..12}.rna.fna.gz \
    https://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.files.installed
wget -q\
    https://github.com/SACGF/cdot/releases/download/v$CDOT_RELEASE/$CDOT_FILENAME
