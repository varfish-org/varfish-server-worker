#!/usr/bin/env bash

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

set -x

mkdir -p $SEQREPO_ROOT_DIR
seqrepo init -i master

# seqrepo load is too verbose and we cannot silence it
2>&1 seqrepo \
    load --instance-name master --namespace ENSEMBL \
    $DATA_DIR/tmp/$GENOME_RELEASE/Homo_sapiens.$ENSEMBL_TOKEN.cdna.all.fa.gz \
| tail

# seqrepo load is too verbose and we cannot silence it
2>&1 seqrepo \
    load --instance-name master --namespace RefSeq \
    $DATA_DIR/tmp/$GENOME_RELEASE/human.{1..12}.rna.fna.gz \
| tail
