# VarFish Server Worker

This repository contains the worker used by VarFish Server to execute certain background task.
They are written in the Rust programming language to speed up the execution of certain tasks.
At the moment, the following worker exist:

- `sv bg-db-to-bin` -- convert background TSV file (from VarFish DB Downloader or `sv inhouse-db-build` to binary format as input for `sv query`)
- `sv inhouse-db-build` -- combine multiple SV TSV files from VarFish Annotator annotated into a background TSV file for the worker
- `sv query` -- perform SV (overlap) annoation and filtration

## Overall Design

For running queries, the worker tool is installed into the VarFish Server image and are run as executables.
The VarFish Server Celery task writes out file(s) with the input for the worker.
The worker then reads in this file, may use some static database files but also the VarFish server postgres database, and can write out out a result file or store the results in the postgres database
The celery worker will then process any worker output further.
The book keeping is done by the Celery worker that runs from the code of VarFish server.

## The `sv bg-db-to-bin` Command

Convert background database from `sv inhouse-db-build` or VarFish DB Downloader output to a binary file.
The binary file uses `flatbuffers` for fast I/O (startup time is less than 200ms on workstation with NVME).

```
$ varfish-server-worker sv bg-db-to-bin \
    --path-input-tsv INPUT.tsv \
    --path-output-bin OUTPUT.bin \
    --input-type TYPE
```

## The `sv inhouse-db-build` Command

Compute in-house database TSV file from one or more VarFish Annotator annotated SV call sets.
You can use `@PATH.txt` for a file that has more TSV files, one path per line.

```
$ varfish-server-worker sv inhouse-db-build \
    --path-output-tsv OUT.tsv \
    INPUT.tsv [INPUT.tsv ...]
```

## The `sv query` Command

This is for actually running the queries.

```
$ varfish-server-worker sv query \
    --path-db PATH_DB \
    --path-query-json QUERY.json \
    --path-input-svs VARFISH_ANNOTATED.gts.tsv.gz
```

## Development Setup

You will need a recent version of flatbuffers, e.g.:

```
# bash utils/install-flatbuffers.sh
# export PATH=$PATH:$HOME/.local/share/flatbuffers/bin
```