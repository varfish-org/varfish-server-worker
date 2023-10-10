#!/usr/bin/bash

# Setup database excerpts required for running the SV queries.
#
# Define `SRC` with the path to the "output/full" directory of the
# `varfish-db-downloader`.
#
# This will pull all data from chromosome 16.

set -x
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
STRUCVARS_AGGREGATE="cargo run --release -- strucvars aggregate"
TXT_TO_BIN="cargo run --release -- strucvars txt-to-bin --assembly=grch37"

mkdir -p $SCRIPT_DIR/db/worker/grch37/{features,genes,strucvars/bgdbs,tads}
mkdir -p $SCRIPT_DIR/db/worker/noref/genes
mkdir -p $SCRIPT_DIR/db/mehari/grch37

# Filter BED files
filter-bed()
{
    zcat --force "$@" \
    | awk -F $'\t' '(($1 ~ /^#/) || ($1 ~ /16$/))'
}

# Write to a file (helpful in interpreting `set -x` output.)
write-to()
{
    tee $1 >/dev/null
}

# mehari/grch37/

mkdir -p /tmp/mehari-data-grch37

pushd $SCRIPT_DIR/mehari-tx-scripts ; \
DATA_DIR=/tmp/mehari-data-grch37 \
ACTIVATE_CONDA="conda activate" \
USE_CONDA_ENV=varfish-server-worker-query-bootstrap \
GENOME_RELEASE=grch37 \
bash run.sh ; \
popd

cp /tmp/mehari-data-grch37/pass-2/txs.bin.zst* \
$SCRIPT_DIR/db/mehari/grch37

# worker/noref/genes/

cp $SRC/mehari/genes-xlink-20230913/genes-xlink.tsv \
$SCRIPT_DIR/db/worker/noref/genes/xlink.tsv
cp $SRC/worker/acmg-sf-3.1+0.10.1/acmg_sf.tsv \
$SCRIPT_DIR/db/worker/noref/genes/acmg.tsv
cp $SRC/worker/mim2gene-20230913+0.10.1/mim2gene.tsv \
$SCRIPT_DIR/db/worker/noref/genes/mim2gene.tsv

# grch37/features/

filter-bed $SRC/tracks/track-features-ucsc-rmsk-grch37-20200322+0/rmsk.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/features/masked_repeat.bed
$TXT_TO_BIN \
    --input-type=masked-region \
    --path-input=$SCRIPT_DIR/db/worker/grch37/features/masked_repeat.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/features/masked_repeat.bin

filter-bed $SRC/tracks/track-features-ucsc-genomicsuperdups-grch37-20111025+0/genomicSuperDups.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/features/masked_seqdup.bed
$TXT_TO_BIN \
    --input-type=masked-region \
    --path-input=$SCRIPT_DIR/db/worker/grch37/features/masked_seqdup.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/features/masked_seqdup.bin

# grch37/strucvars/

filter-bed tests/strucvars/query/Case_3_ingested.vcf \
| $STRUCVARS_AGGREGATE --path-output $SCRIPT_DIR/db/worker/grch37/strucvars/inhouse.bed /dev/stdin
$TXT_TO_BIN \
    --input-type=strucvar-inhouse \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/inhouse.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/inhouse.bin

$TXT_TO_BIN \
    --input-type=clinvar-sv \
    --path-input=$SCRIPT_DIR/db/worker/grch37/clinvar-full-release.head.jsonl.gz \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/clinvar.bin

filter-bed $SRC/tracks/track-strucvars-dbvar-grch37-20230516+0/dbvar.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dbvar.bed
$TXT_TO_BIN \
    --input-type=strucvar-db-var \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dbvar.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dbvar.bin

filter-bed $SRC/tracks/track-strucvars-dgv-grch37-20200225+0/dgv.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dgv.bed
$TXT_TO_BIN \
    --input-type=strucvar-dgv \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dgv.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dgv.bin

filter-bed $SRC/tracks/track-strucvars-dgv-gs-grch37-20200225+0/dgv-gs.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dgv-gs.bed
$TXT_TO_BIN \
    --input-type=strucvar-dgv-gs \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dgv-gs.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/dgv-gs.bin

filter-bed $SRC/tracks/track-strucvars-exac-grch37-0.3.1+0/exac.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/exac.bed
$TXT_TO_BIN \
    --input-type=strucvar-exac-cnv \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/exac.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/exac.bin

filter-bed $SRC/tracks/track-strucvars-g1k-grch37-phase3v2+0/g1k.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/g1k.bed
$TXT_TO_BIN \
    --input-type=strucvar-g1k \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/g1k.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/g1k.bin

filter-bed $SRC/tracks/track-strucvars-gnomad-grch37-2.1.1+0/gnomad.bed.gz \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/gnomad.bed
$TXT_TO_BIN \
    --input-type=strucvar-gnomad-sv \
    --path-input=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/gnomad.bed \
    --path-output=$SCRIPT_DIR/db/worker/grch37/strucvars/bgdbs/gnomad.bin

# grch37/tads/

filter-bed $SRC/worker/tads-grch37-dixon2015/hesc.bed \
| write-to $SCRIPT_DIR/db/worker/grch37/tads/hesc.bed

filter-bed $SRC/worker/patho-mms-grch37-20220730+0.10.1/patho-mms.bed \
| write-to $SCRIPT_DIR/db/worker/grch37/strucvars/patho-mms.bed
