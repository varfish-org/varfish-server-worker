#!/usr/bin/env bash

# Build mehari transcript protobuf file (pass 1).
#
# Because the NCBI transcripts are not stable, fetching the transcript sequences
# will fail for certain transcripts.  This will be apparent from the report file
# written out by this pass.  The `fix-seqrepo.sh` script will attempt to work
# around this issue by fetching the sequences directly from NCBI and then `pass-2.sh`
# will work without such failures.

set -euo pipefail

export USE_CONDA_ENV=${USE_CONDA_ENV-mehari-data}
export DATA_DIR=${DATA_DIR-/tmp/mehari-data}
export SEQREPO_ROOT_DIR=$DATA_DIR/seqrepo
export PATH=$PATH:/home/runner/micromamba-root/envs/mehari-data/bin

set -x

mkdir -p $DATA_DIR/pass-1

symbols="$(echo ACSF3 ACSM2B ACSM3 ADHD1 ARL6IP1 ARMC5 BMIQ5 C16orf58 C16orf71 C16orf82 C16orf84 C16orf95 C16orf96 CARHSP1 CASP16P CCDC113 Ccdc78 CDIPT CFDP1 CHDS1 CIAPIN1 CKLF CLUAP1 CMTM2 CCDC135 COTL1 CPNE7 CTRL DCTPP1 DHX38 EMP2 ENKD1 ERAF FAHD1 FAM57B FBRS FOXC2-AS1 GLG1 HBAP1 HIRIP3 HN1L IBD8 IHPS2 ITFG3 KDM8 KIAA0895L LOC124220 LOC81691 LUC7L LYPLA3 MC1R MCOPCT1 METRN MKL2 MPHOSPH6 MT1G MT1X NIP30 NOB1 NOMO1 NPW NUBP2 NUPR1 OGFOD1 PDF PDPR PKDTS PMFBP1 POLR3K PRMT7 PRR35 RPS15A RSL1D1 SHCBP1 SLZ1 SNAI3-AS1 SNORD71 SPSB3 SRCAP TANGO6 TAO2 TBC1D24 TEDC2 TELO2 TMEM112 TMEM8A TNRC6A TSR3 UNKL VAT1L VPS35L WFDC1 ZG16 ZNF23 ZNF200 ZNF263 ZNF629 ZNF843 | tr ' ' '\n')"

mehari db create txs \
    --path-out $DATA_DIR/pass-1/txs.bin.zst \
    --path-seqrepo-instance $DATA_DIR/seqrepo/master \
    --path-cdot-json $DATA_DIR/tmp/$GENOME_RELEASE/$CDOT_FILENAME \
    --genome-release $GENOME_RELEASE \
    $(for symbol in $symbols; do echo --gene-symbols=$symbol; done)

# Ensure that the output can be decompressed.
zstd -c -d $DATA_DIR/pass-1/txs.bin.zst > /dev/null
