[![Crates.io](https://img.shields.io/crates/d/varfish-server-worker.svg)](https://crates.io/crates/varfish-server-worker)
[![Crates.io](https://img.shields.io/crates/v/varfish-server-worker.svg)](https://crates.io/crates/varfish-server-worker)
[![Crates.io](https://img.shields.io/crates/l/varfish-server-worker.svg)](https://crates.io/crates/varfish-server-worker)
[![CI](https://github.com/bihealth/varfish-server-worker/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/varfish-server-worker/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/varfish-server-worker/branch/main/graph/badge.svg?token=t5oheMdukf)](https://codecov.io/gh/bihealth/varfish-server-worker)
[![DOI](https://zenodo.org/badge/590461457.svg)](https://zenodo.org/badge/latestdoi/590461457)

# VarFish Server Worker

This repository contains the worker used by VarFish Server to execute certain background task.
They are written in the Rust programming language to speed up the execution of certain tasks.
At the moment, the following sub commands exist:

- `db compile` -- compile files from `varfish-db-downloader` to be useable by the query sub commands
- `db mk-inhouse` -- compile per-case structural variant into an in-house database previously created by `db compile`
- `sv query` -- perform SV (overlap) annoation and filtration

## Overall Design

For running queries, the worker tool is installed into the VarFish Server image and are run as executables.
The VarFish Server Celery task writes out file(s) with the input for the worker.
The worker then reads in this file, may use some static database files but also the VarFish server postgres database, and can write out out a result file or store the results in the postgres database.
The celery worker will then process any worker output further.
The book keeping is done by the Celery worker that runs from the code of VarFish server.

Future versions may provide persistently running HTTP/REST servers that provide functionality without startup cost.

## The `db compile` Command

Convert output of `varfish-db-downloader` to a directory with databases to be used by query commands such as `sv query`.

```
$ varfish-server-worker db compile \
    --path-db-downloader SRC/varfish-db-downloader \
    --path-worker-db DST/varfish-server-worker-db
```

## The `db mk-inhouse` Command

Import multiple files created by `varfish-annotator annotate_svs` into a database previously created by `db compile`.
You can specify the files individually.
Paths starting with an at (`@`) character are interpreted as files with lists of paths.
You can mix paths with `@` and without.

```
$ varfish-server-worker db mk-inhouse \
    --path-worker-db WORK/varfish-server-worker-db \
    IN/file1.gts.tsv.gz [IN/file1.gts.tsv.gz] \

# OR:

$ varfish-server-worker db mk-inhouse \
    --path-worker-db WORK/varfish-server-worker-db \
    @IN/path-list.txt [@IN/path-list2.txt]
```

## The `sv query` Command

This is for actually running the queries.

```
$ varfish-server-worker sv query \
    --path-db PATH_DB \
    --path-query-json QUERY.json \
    --path-input-svs VARFISH_ANNOTATED.gts.tsv.gz
```

## The `server rest` Command

This runs the REST API.


```
$ varfish-server-worker server rest \
    --path-db PATH_DB
```

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
$ terraform import github_repository.reev reev
$ terraform validate
$ terraform fmt
$ terraform plan
$ terraform apply
```
