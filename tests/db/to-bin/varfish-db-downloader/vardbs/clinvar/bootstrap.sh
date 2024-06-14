#!/usr/bin/bash

RELEASE_DATE=${RELEASE_DATE-20240612}
CLINVAR_THIS_VERSION=${CLINVAR_THIS_VERSION-0.17.0}
VERSION=${RELEASE_DATE}+${CLINVAR_THIS_VERSION}

URL=https://github.com/varfish-org/clinvar-data-jsonl/releases/download/clinvar-weekly-${RELEASE_DATE}/clinvar-data-extract-vars-${VERSION}.tar.gz

RELEASE=grch37

# # Download subset from clinvar-data-jsonl for subsetting.

set -x
set -euo pipefail

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

export TMPDIR=$(mktemp -d)
#trap "rm -rf $TMPDIR" EXIT

DST=$TMPDIR/clinvar-data-extract-vars-${VERSION}.tar.gz
wget -O $DST $URL
tar -C $TMPDIR -xf $DST

zcat $TMPDIR/clinvar-data-extract-vars-${VERSION}/clinvar-variants-${RELEASE}-strucvars.jsonl.gz \
| head -n 500 \
| tee $SCRIPT_DIR/clinvar-svs.jsonl \
| gzip -c \
> $SCRIPT_DIR/clinvar-svs.jsonl.gz
