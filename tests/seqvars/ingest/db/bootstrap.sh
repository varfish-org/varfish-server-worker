#!/usr/bin/bash

# Helper script that bootstraps the database for the tests.
#
# We include transcripts, clinvar, and frequencies for BRCA1 only, clinvar
# and frequencies are augmented with chrMT.

set -x

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

export TMPDIR=$(mktemp -d)

# Select BRCA1 and chrMT
echo -e "17\t41176014\t41297797" > $TMPDIR/brca1.grch37.bed
echo -e "chr17\t43024028\t43145632" > $TMPDIR/brca1.grch38.bed
echo -e "MT\t0\t16569" >> $TMPDIR/brca1.grch37.bed
echo -e "chrMT\t0\t16569" >> $TMPDIR/brca1.grch37.bed

CLINVAR_DATE=20240612
CLINVAR_ANNONARS=0.38.0
MEHARI_TXS=0.4.4
S5CMD="s5cmd --endpoint-url=https://ceph-s3-public.cubi.bihealth.org --no-sign-request --no-verify-ssl"
MEHARI_FREQS_VERSION_37=2.1.1+2.1.1+3.1+20200327+0.33.0
MEHARI_FREQS_VERSION_38=4.0+4.0+3.1+20200327+0.33.0
GENES="BRCA1 OPA1 AKR1C3 TTN"

for release in grch37 grch38; do
    # mehari clinvar
    if [[ -e $SCRIPT_DIR/$release/seqvars/clinvar/rocksdb ]]; then
        >&2 echo "WARNING: $SCRIPT_DIR/$release/seqvars/clinvar already exists"
    else
        mkdir -p $SCRIPT_DIR/$release/seqvars/clinvar
        wget -O $TMPDIR/annonars-clinvar-minimal-$release-$CLINVAR_DATE+$CLINVAR_ANNONARS.tar.gz \
            https://github.com/varfish-org/annonars-data-clinvar/releases/download/annonars-data-clinvar-$CLINVAR_DATE/annonars-clinvar-minimal-$release-$CLINVAR_DATE+$CLINVAR_ANNONARS.tar.gz
        tar -C $TMPDIR -xf $TMPDIR/annonars-clinvar-minimal-$release-$CLINVAR_DATE+$CLINVAR_ANNONARS.tar.gz

        annonars db-utils copy \
            --skip-cfs clinvar_by_accession \
            --path-in $TMPDIR/annonars-clinvar-minimal-$release-$CLINVAR_DATE+$CLINVAR_ANNONARS/rocksdb \
            --path-out $SCRIPT_DIR/$release/seqvars/clinvar/rocksdb \
            --path-beds $TMPDIR/brca1.$release.bed
    fi

    # mehari freqs
    if [[ -e $SCRIPT_DIR/$release/seqvars/freqs/rocksdb ]]; then
        >&2 echo "WARNING: $SCRIPT_DIR/$release/seqvars/freqs already exists"
    else
        mkdir -p $SCRIPT_DIR/$release/seqvars/freqs/rocksdb
        if [[ $release == "grch37" ]]; then
            $S5CMD sync \
                "s3://varfish-public/reduced-dev/mehari/freqs-grch37-${MEHARI_FREQS_VERSION_37}/rocksdb/*" \
                $SCRIPT_DIR/$release/seqvars/freqs/rocksdb
        else
            $S5CMD sync \
                "s3://varfish-public/reduced-dev/mehari/freqs-grch38-${MEHARI_FREQS_VERSION_38}/rocksdb/*" \
                $SCRIPT_DIR/$release/seqvars/freqs/rocksdb
        fi
    fi

    # tx database
    if [[ -e $SCRIPT_DIR/$release/txs.bin.zst ]]; then
        >&2 echo "WARNING: $SCRIPT_DIR/$release/txs.bin.zst already exists"
    else
        mkdir -p $SCRIPT_DIR/$release
        wget -O $TMPDIR/mehari-data-txs-$release-$MEHARI_TXS.bin.zst \
            https://github.com/bihealth/mehari-data-tx/releases/download/v$MEHARI_TXS/mehari-data-txs-$release-$MEHARI_TXS.bin.zst
        mehari db subset \
            --path-in $TMPDIR/mehari-data-txs-$release-$MEHARI_TXS.bin.zst \
            --path-out $SCRIPT_DIR/$release/txs.bin.zst \
            $(for gene in $GENES; do echo --gene-symbols $gene; done)
    fi
done
