[![Crates.io](https://img.shields.io/crates/d/varfish-server-worker.svg)](https://crates.io/crates/varfish-server-worker)
[![CI](https://github.com/bihealth/varfish-server-worker/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/varfish-server-worker/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/varfish-server-worker/branch/main/graph/badge.svg?token=t5oheMdukf)](https://codecov.io/gh/bihealth/varfish-server-worker)
[![DOI](https://zenodo.org/badge/590461457.svg)](https://zenodo.org/badge/latestdoi/590461457)

# VarFish Server Worker

> [!NOTE]
> This repository contains code that runs inside a VarFish Server.
> If you are looking to run your own VarFish Server, look [here at bihealth/varfish-server](https://github.com/bihealth/varfish-server).

This repository contains the worker used by VarFish Server to execute certain background task.
They are written in the Rust programming language to speed up the execution of certain tasks.
At the moment, the following sub commands exist:

- `db` -- subcommands to build binary (protobuf) database files
- `seqvars` -- subcommands for processing sequence (aka small/SNV/indel) variants
    - `seqvars ingest` -- convert single VCF file into internal format for use with `seqvars query`
    - `seqvars query` -- perform sequence variant filtration and on-the-fly annotation
    - `seqvars prefilter` -- limit the result of `seqvars prefilter` by population frequency and/or distance to exon
    - `seqvars aggregate` -- read through multiple VCF files written by `seqvars ingest` and computes a carrier counts table.
- `strucvars` -- subcommands for processing structural (aka large variants, CNVs, etc.) variants
    - `strucvars ingest` -- convert one or more structural variant files for use with `strucvars query`
    - `strucvars aggregate` -- compile per-case structural variant into an in-house database, to be converted to `.bin` with `strucvars txt-to-bin`.
    - `strucvars txt-to-bin` -- convert text files downloaded by [varfish-db-downloader](https://github.com/bihealth/varfish-db-downloader/) to binary for fast use in `strucvars query` commands
    - `strucvars query` -- perform structural variant filtration and on-the-fly annotation

## Overall Design

For running queries, the worker tool is installed into the VarFish Server image and are run as executables.
Internally, VarFish server works on VCF files stored in an S3 storage.

For import, the user gives the server access to the VCF files to import.
The server will then use the worker executable to ingest the data into the internal format using `{seqvars,strucvars} ingest`.
These files are then stored in the internal S3 storage.

For queries, the server will create a query JSON file and then pass this query JSON file together with the internal file to the worker executable.
The worker will create a result file that can be directly imported by the server to be displayed to the user.

Future versions may provide persistently running HTTP/REST servers that provide functionality without startup cost.

## The `seqvars ingest` Command

This command takes as the input a single VCF file from a (supported) variant caller and converts it into a file for further querying.
The command interprets the following fields which are written out by the commonly used variant callers such as GATK UnifiedGenotyper, GATK HaplotypeCaller, and Illumina Dragen.

- `FORMAT/GT` -- genotype
    - the following `GT` values are written out as `0/0`, `0/1`, `1/0`, `1/1`, `0|0`, `0|1`, `1|0`, `1|1`, `./.`, `.|.`, `.`
    - no combination of no-call (`.`) and called allele is written out
- `FORMAT/GQ` -- genotype quality
- `FORMAT/DP` -- total read coverage
- `FORMAT/AD` -- allelic depth, one value per allele (including reference0)
- `FORMAT/PS` -- physical phasing information as written out by GATK HaplotypeCaller in GVCF workflow and Dragen variant caller
- `FORMAT/SQ` -- "somatic quality" for each alternate allele, as written out by Illumina Dragen variant caller
    - this field will be written as `FORMAT/GQ`

The `seqvars ingest` command will annotate the variants with the following information:

- gnomAD genomes and exomes allele frequencies
- gnomAD-mtDNA and HelixMtDb allele frequencies
- functional annotation following the [VCF ANN field standard](https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf)
    - `Gene_Name` is writen as HGNC symbol
    - `Gene_ID` is written as HGNC ID

The command will emit one output line for each variant allele from the input and each affected gene.
That is, if two variant alleles affect two genes, four records will be written to the output file.
The annotation will be written out for one highest impact.

Overall, the command will emit the following header rows in addition to the `##contig=<ID=.,length=.>` lines.

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=gnomad_exomes_an,Number=1,Type=Integer,Description="Number of alleles in gnomAD exomes">
##INFO=<ID=gnomad_exomes_hom,Number=1,Type=Integer,Description="Number of hom. alt. carriers in gnomAD exomes">
##INFO=<ID=gnomad_exomes_het,Number=1,Type=Integer,Description="Number of het. alt. carriers in gnomAD exomes">
##INFO=<ID=gnomad_exomes_hemi,Number=1,Type=Integer,Description="Number of hemi. alt. carriers in gnomAD exomes">
##INFO=<ID=gnomad_genomes_an,Number=1,Type=Integer,Description="Number of alleles in gnomAD genomes">
##INFO=<ID=gnomad_genomes_hom,Number=1,Type=Integer,Description="Number of hom. alt. carriers in gnomAD genomes">
##INFO=<ID=gnomad_genomes_het,Number=1,Type=Integer,Description="Number of het. alt. carriers in gnomAD genomes">
##INFO=<ID=gnomad_genomes_hemi,Number=1,Type=Integer,Description="Number of hemi. alt. carriers in gnomAD genomes">
##INFO=<ID=helix_an,Number=1,Type=Integer,Description="Number of alleles in HelixMtDb">
##INFO=<ID=helix_hom,Number=1,Type=Integer,Description="Number of hom. alt. carriers in HelixMtDb">
##INFO=<ID=helix_het,Number=1,Type=Integer,Description="Number of het. alt. carriers in HelixMtDb">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##x-varfish-case-uuid=d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c
##x-varfish-version=<ID=varfish-server-worker,Version=x.y.z>
##x-varfish-version=<ID=orig-caller,Name=Dragen,Version=SW: 07.021.624.3.10.9, HW: 07.021.624>
##x-varfish-version=<ID=orig-caller,Name=GatkHaplotypeCaller,Version=4.4.0.0>
```

> [!NOTE]
> The gnomad-mtDNA information is written to the `INFO/gnomdad_genome_*` fields.

> [!NOTE]
> Future versions of the worker will annotate the worst effect on a MANE select or MANE Clinical transcript.

## The `seqvars prefilter` Command

This file takes as the input a file created by `seqvars ingest` and filters the variants by population frequency and/or distance to exon.
You can pass the prefilter criteria as JSON on the command line corresponding to the following Rust structs:

```rust
struct PrefilterParams {
    /// Path to output file.
    pub path_out: String,
    /// Maximal allele population frequency.
    pub max_freq: f64,
    /// Maximal distance to exon.
    pub max_dist: i32,
}
```

You can either specify the parameters on the command line directly or pass a path to a JSONL file starting with `@`.
You can mix both ways.

```
$ varfish-server-worker strucvars prefilter \
    --path-input INPUT.vcf \
    --params '{"path_out": "out.vcf", "max_freq": 0.01, "max_dist": 100}' \
    [--params ...] \

# OR

$ varfish-server-worker strucvars prefilter \
    --path-input INPUT.vcf \
    --params @path/to/params.json \
    [--params ...] \

## The `seqvars aggregate` Command

This command reads through multiple files written by `seqvars ingest` and computes a in-house carrier counts table.
You can specify the VCF files individually on the command line or pass in files that have paths to the VCF files line by line.
The resulting table is a folder to a RocksDB database.

```shell session
varfish-server-worker seqvars aggregate \
    --genome-build {grch37,grch38} \
    --path-out-rocksdb rocksdb/folder \
    --path-in-vcf path/to/vcf.gz \
    --path-in-vcf @path/to/file/list.txt
```

## The `seqvars query` Command

This command perform the querying of sequence variants and further annotation using annonars databases.

## The `strucvars ingest` Command

This command takes as the input one or more VCF files from structural variant callers and converts it into a file for further querying.
The command supports the following variant callers and can guess the caller from the VCF header and first record.

- Delly2
- Dragen-SV (equivalent to Manta)
- Dragen-CNV
- GATK gCNV
- Manta
- MELT
- PopDel
- Sniffles2

One record will be written out for each variant, each with a single alternate allele.

The following symbolic `ALT` alleles are used:

- `<DEL>`
- `<DUP>`
- `<INS>`
- `<INV>`
- VCF break-end syntax, e.g., `T[chr1:5[`

The following `INFO` fields are written:

- `IMPRECISE` -- flag that specifies that this is an imprecise variant
- `END` -- end position of the variants
- `SVTYPE` -- type of the variant, one of `<DEL>`, `<DUP>`, `<INS>`, `<INV>`, `BND`
- `SVLEN` -- absolute length of the SV for linear variants, `.` for non-linear variants
- `SVCLAIM` -- specificaton of `D` (change in abundance), `J` (novel junction), or `DJ` (both change in abundance and novel junction)
- `callers` -- (non-standard field), list of callers that called the variant
- `chr2` -- (non-standard field), second chromosome for BND variants
- `annsv` -- (non-standard field), annotation of the variant effect on each affected gene

The `annsv` field is a pipe-character (`|`) separated list of the following fields:

1. symbolic alternate alele, e.g., `<DEL>`
2. effects on the gene's transcript, separated by `&`
    - `transcript_variant` -- variant affects the whole transcript
    - `exon_variant` -- variant affects exon
    - `splice_region_variant` -- variant affects splice region
    - `intron_variant` -- variant affects only intron
    - `upstream_variant` -- variant upsream of gene
    - `downstream_variant` -- variant downstream of gene
    - `intergenic_variant` -- default for "no gene affected", but never written
3. HGNC gene symbol, e.g., `BRCA1`
4. HGNC gene ID, e.g., `HGNC:1100`

The following `FORMAT` fields are written:

- `GT` -- (standard field) genotype, if applicable
- `GQ` -- (standard field) genotype quality, if applicable
- `pec` -- total coverage with paired-end reads
- `pev` -- paired-end reads supporting the variant
- `src` -- total coverage with split reads
- `srv` -- split reads supporting the variant
- `amq` -- average mapping quality over the variant
- `cn` -- copy number of the variant in the sample
- `anc` -- average normalized coverage over the variant in the sample
- `pc` -- point count (windows/targets/probes)


Overall, the command will emit the following header rows in addition to the `##contig=<ID=.,length=.>` lines.

```
##fileformat=VCFv4.4
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Length of structural variant">
##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively">
##INFO=<ID=callers,Number=.,Type=String,Description="Callers that called the variant">
##INFO=<ID=chr2,Number=1,Type=String,Description="Second chromosome, if not equal to CHROM">
##INFO=<ID=annsv,Number=1,Type=String,Description="Effect annotations: 'Allele | Annotation | Gene_Name | Gene_ID'">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=pec,Number=1,Type=Integer,Description="Total coverage with paired-end reads">
##FORMAT=<ID=pev,Number=1,Type=Integer,Description="Paired-end reads supporting the variant">
##FORMAT=<ID=src,Number=1,Type=Integer,Description="Total coverage with split reads">
##FORMAT=<ID=srv,Number=1,Type=Integer,Description="Split reads supporting the variant">
##FORMAT=<ID=amq,Number=1,Type=Float,Description="Average mapping quality over the variant">
##FORMAT=<ID=cn,Number=1,Type=Integer,Description="Copy number of the variant in the sample">
##FORMAT=<ID=anc,Number=1,Type=Float,Description="Average normalized coverage over the variant in the sample">
##FORMAT=<ID=pc,Number=1,Type=Integer,Description="Point count (windows/targets/probes)">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=CNV,Description="Copy Number Variation">
##ALT=<ID=INV,Description="Inversion">
##fileDate=20230421
##x-varfish-genome-build=GRCh37
##SAMPLE=<ID=index,Sex="Male",Disease="Affected">
##SAMPLE=<ID=father,Sex="Male",Disease="Unaffected">
##SAMPLE=<ID=mother,Sex="Female",Disease="Unaffected">
##PEDIGREE=<ID=index,Father="father",Mother="mother">
##PEDIGREE=<ID=father>
##PEDIGREE=<ID=mother>
##x-varfish-case-uuid=d2bad2ec-a75d-44b9-bd0a-83a3f1331b7c
##x-varfish-version=<ID=varfish-server-worker,Version="x.y.z">
##x-varfish-version=<ID=Delly,Name="Delly",Version="1.1.3">
```

> [!NOTE]
> The `strucvars ingest` step does not perform any annotation.
> It only merges the input VCF files from multiple callers (all files must have the same samples) and converts them into the internal format.
> The `INFO/annsv` field is filled by `strucvars query`.

## The `strucvars aggregate` Command

Import multiple files created by `strucvars ingest` into a database that can be convered to `.bin` with `strucvars txt-to-bin` and then used by `strucvars query`.
You can specify the files individually.
Paths starting with an at (`@`) character are interpreted as files with lists of paths.
You can mix paths with `@` and without.

```
$ varfish-server-worker strucvars aggregate \
    --genome-release {Grch37,Grch38} \
    --path-output OUT.tsv \
    --path-input IN/file1.vcf.gz \
    [--path-input IN/file1.vcf.gz] \

# OR:

$ varfish-server-worker db mk-inhouse \
    --genome-release {Grch37,Grch38} \
    --path-output OUT.tsv \
    --path-input @IN/path-list.txt \
    [--path-input @IN/path-list2.txt]
```

## The `strucvars txt-to-bin` Command

Convert output of [varfish-db-downloader](https://github.com/bihealth/varfish-db-downloader/) to a directory with databases to be used by query commands such as `strucvars query`.

```
$ varfish-server-worker strucvars txt-to-bin \
    --input-type {ClinvarSv,StrucvarInhouse,...} \
    --path-input IN.txt \
    --path-output DST.bin
```

## The `strucvars query` Command

Run a query on a VCF file with structural variants as created by `strucvars ingest` using a varfish worker database.

```
$ varfish-server-worker strucvars query \
    --genome-release grch37 \
    --path-db path/to/worker-db \
    --path-input IN.vcf.gz \
    --path-output OUT.jsonl
```

The worker database has the following structure:

```
$ROOT/
    noref/
        genes/
            acmg.tsv              -- ACMG SF list genes
            mim2gene.tsv          -- OMIM to NCBI mapping from clingen
            xlink.bin             -- gene crosslinks
    {genome_release}/             -- one per genome release
        mehari/
            txs.bin.zstd          -- mehari transcripts
        features/                 -- features important for annotation
            masked_repeat.bin     -- masked repeats
            masked_seqdup.bin     -- masked segmental duplications
        strucvars/                -- structural variant specific
            bgdbs/                -- background databases
                dbvar.bin         -- dbVar
                dgv.bin           -- DGV
                dgv-gs.bin        -- DGV gold standard
                exac.bin          -- ExAC CNVs
                g1k.bin           -- 1000 genomes CNVs
                gnomad.bin        -- gnomAD-SVs
            clinvar.bin           -- ClinVar SVs
            inhouse.bin           -- inhouse SV database
            patho-mms.bed         -- well-known pathogenic DELs/DUPs
        tads/
            hesc.bed              -- hESC TAD definitions
```

# Developer Information

This section is only relevant for developers of `varfish-server-worker`.

## Development Setup

You will need a recent version of protocolbuffers, e.g.:

```
# bash utils/install-protoc.sh
# export PATH=$PATH:$HOME/.local/share/protoc/bin
```

## GitHub Project Management

We use Terraform for managing the GitHub project settings (as applicable):

```
$ export GITHUB_OWNER=bihealth
$ export GITHUB_TOKEN=ghp_<thetoken>

$ cd utils/terraform
$ terraform init
$ terraform import github_repository.varfish-sever-worker varfish-sever-worker
$ terraform validate
$ terraform fmt
$ terraform plan
$ terraform apply
```
