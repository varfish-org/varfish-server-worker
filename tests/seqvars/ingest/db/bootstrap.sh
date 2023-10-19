#!/usr/bin/bash

# Helper script that bootstraps the database for the tests.
#
# We include transcripts, clinvar, and frequencies for BRCA1 only.

set -x

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p $SCRIPT_DIR/{grch37,grch38}

# transcripts

wget -O /tmp/cdot-0.2.12.refseq.grch37_grch38.json.gz \
    https://github.com/SACGF/cdot/releases/download/v0.2.12/cdot-0.2.12.refseq.grch37_grch38.json.gz

mehari \
    db \
    create \
    txs \
    --path-out \
    $SCRIPT_DIR/grch37/txs.bin.zst \
    --path-cdot-json \
    /tmp/cdot-0.2.12.refseq.grch37_grch38.json.gz \
    --path-seqrepo-instance \
    ~/Data/mehari/db-workspace/b37/seqrepo-data/master \
    --genome-release \
    grch37 \
    --gene-symbols BRCA1

# clinvar

wget -O /tmp/annonars-clinvar-minimal-grch37-20231015+0.24.1.tar.gz \
    https://github.com/bihealth/annonars-data-clinvar/releases/download/annonars-data-clinvar-20231015/annonars-clinvar-minimal-grch37-20231015+0.24.1.tar.gz
tar -C /tmp -xvf /tmp/annonars-clinvar-minimal-grch37-20231015+0.24.1.tar.gz

mkdir -p $SCRIPT_DIR/grch37/seqvars/clinvar
cp /tmp/annonars-clinvar-minimal-grch37-20231015+0.24.1/spec.yaml \
    $SCRIPT_DIR/grch37/seqvars/clinvar
annonars db-utils copy \
    --path-in /tmp/annonars-clinvar-minimal-grch37-20231015+0.24.1/rocksdb \
    --path-out $SCRIPT_DIR/grch37/seqvars/clinvar/rocksdb \
    --range 17:41183866:41337086

# frequencies

DOWNLOADER=/data/sshfs/data/cephfs-1/work/projects/cubi_varfish_data/2023-06-05_varfish-db-downloader/varfish-db-downloader
mkdir -p mkdir -p $SCRIPT_DIR/grch37/seqvars/freqs
cp $DOWNLOADER/output/full/mehari/freqs-grch37-2.1.1+2.1.1+3.1+20200327+0.19.0/spec.yaml \
    $SCRIPT_DIR/grch37/seqvars/freqs
annonars db-utils copy \
    --path-in $DOWNLOADER/output/full/mehari/freqs-grch37-2.1.1+2.1.1+3.1+20200327+0.19.0/rocksdb \
    --path-out $SCRIPT_DIR/grch37/seqvars/freqs/rocksdb \
    --range 17:41183866:41337086
