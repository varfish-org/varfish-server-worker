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
# sv-query --result-set-id INPUT.vcf.gz OUTPUT.tsv
```

**Param: `--result-set-id` (required)**

This is written out to the column `svqueryresultset`.

**Input: `INPUT.vcf.gz`**

A VCF file with the structural variants.
This is dumped from the 

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