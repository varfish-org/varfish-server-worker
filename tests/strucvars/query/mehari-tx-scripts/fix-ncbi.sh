#!/usr/bin/env bash

# cf. http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

set -x

mkdir -p $DATA_DIR/for-fix

{ grep "because it has no sequence" $DATA_DIR/pass-1/txs.bin.zst.report || true; } \
| cut -f 2 \
> $DATA_DIR/for-fix/missing.txt

wc -l $DATA_DIR/for-fix/missing.txt

split -l 20 $DATA_DIR/for-fix/missing.txt $DATA_DIR/for-fix/missing.txt_

if [[ ! -e $(echo $DATA_DIR/for-fix/missing.txt_* | tr ' ' '\n' | head -n 1) ]]; then
    # nothing to do
    touch $DATA_DIR/for-fix/missing.fasta
    exit 0
fi

for f in $DATA_DIR/for-fix/missing.txt_*; do
    esearch -db nucleotide -query "$(cat $f | tr '\n' ' ')" \
    | efetch -format fasta \
    >> $DATA_DIR/for-fix/missing.fasta
done

2>&1 seqrepo load --instance-name master --namespace RefSeq \
    $DATA_DIR/for-fix/missing.fasta \
| tail
