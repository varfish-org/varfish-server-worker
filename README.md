# VarFish Server Workers

This repository contains the workers used by VarFish Server to execute certain background task.
They are written in the Rust programming language to speed up the execution of certain tasks.
At the moment, the following workers exist:

- `sv-query` -- execution of SV queries

## Overall Design

The workers are installed into the VarFish Server image and are run as executables.
The VarFish Server Celery task writes out file(s) with the input for the worker.
The worker then reads in this file, may use some static database files but also the VarFish server postgres database, and can write out out a result file or store the results in the postgres database
The celery worker will then process any worker output further.
The book keeping is done by the Celery worker that runs from the code of VarFish server.

## The `sv-query` Worker

This worker provides the query functionality.

**Call**

```bash
$ sv-query --result-set-id ID INPUT.vcf.gz OUTPUT.tsv
```

**Param: `--result-set-id ID` (required)**

This is written out to the column `svqueryresultset`.

**Input: `INPUT.vcf.gz`**

A VCF file with the structural variants.
This is dumped from the database by the Celery worker.

The file will have a `##reference=GRCh37` or `##reference=GRCh38` header depending on the reference of the case.

The following `INFO` fields will be written out.

- `INFO/END` -- Integer with end position, if applicable.
- `INFO/POS2` -- Integer with position on second chromosome, if applicable.
- `INFO/CHR2` -- String with name of second chromosome.
- `INFO/CALLER` -- String with caller(s) that found this variant.
- `INFO/SV_TYPE` -- String with SV type, e.g., `DEL`
- `INFO/SV_SUB_TYPE` -- String with SV sub type, e.g., `DEL:ME:SVA`.

The following fields will be writte to the per-sample `FORMAT` fields.

- `FORMAT/GT` -- String with genotype.
- `FORMAT/FT` -- Per-genotype filter values.
- `FORMAT/GQ` -- Float or integer of Phred-scaled genotype quality.
- `FORMAT/PEC` -- Integer, total number of covering paired-end reads.
- `FORMAT/PEV` -- Integer, number of paired-end variant reads.
- `FORMAT/SRC` -- Integer, total number of covering split reads reads.
- `FORMAT/SRV` -- Integer, number of variant split-reads.
- `FORMAT/CN` -- Integer with copy number estimate.
- `FORMAT/ANC` -- Float with average normalized coverage.
- `FORMAT/PC` -- Integer with number of supporting points (for CNVs).

**Output: `OUTPUT.tsv`**
A text file with the query result.
This is a TSV file that can directly be imported into the PostgreSQL database.
The first line contains a header.
Each of the following line contains one result record.

The fields are as follows:

- `sodar_uuid` -- String with automatically generated UUID4 for the result record.
- `svqueryresultset` -- Integer with ID from argument `--result-set-id`.
- `release` -- String with genome build.
- `chromosome` -- String with chromosome name.
- `chromosome_no` -- Integer with number of `chromosome`.
- `bin` -- Integer with UCSC bin of interval or `chromosome` position.
- `chromosome2` -- String with chromosome name, e.g., `X` for GRCh37 and `chrX` for GRCh38.
- `chromosome_no2` -- Integer with chromosome number of `chromsome2`.
- `bin2` -- Integer with UCSC bin of position on `chromosome2`.
- `start` -- Integer with 1-based start position.
- `end` -- Integer with 1-based end position.
- `pe_orientation` -- String with PE orientation, empty or one of `5to3`, `3to5`, `3to3`, or `5to5`.
- `sv_type` -- String with the SV type (e.g., `DEL`)
- `sv_sub_type` -- String with the SV sub type (e.g., `DEL:ME:ALU`)
- `payload` -- JSON field with the result record payload.