#!/usr/bin/bash

mkdir -p features/grch37/{gene_regions,tads}

touch features/grch37/gene_regions/{ensembl,refseq}.spec.json
for path in features/grch37/gene_regions/{ensembl,refseq}{.bed.gz,.bed.gz.tbi}; do
    touch $path
    md5sum $path >$path.md5
done

touch features/grch37/tads/{hesc,imr90}.spec.json
for path in features/grch37/tads/{hesc,imr90}.bed; do
    touch $path
    md5sum $path >$path.md5
done


mkdir -p genes/{acmg,gnomad_constraints,mim2gene,xlink}

touch genes/acmg/acmg.spec.json
touch genes/acmg/acmg.tsv
md5sum genes/acmg/acmg.tsv >genes/acmg/acmg.tsv.md5

touch genes/gnomad_constraints/gnomad_constraints.spec.json
touch genes/gnomad_constraints/gnomad_constraints.tsv
md5sum genes/gnomad_constraints/gnomad_constraints.tsv >genes/gnomad_constraints/gnomad_constraints.tsv.md5

touch genes/mim2gene/mim2gene.spec.json
touch genes/mim2gene/mim2gene.tsv
md5sum genes/mim2gene/mim2gene.tsv >genes/mim2gene/mim2gene.tsv.md5

touch genes/xlink/{ensembl,hgnc}.spec.json
for path in genes/xlink/{ensembl,hgnc}.tsv; do
    touch $path
    md5sum $path >$path.md5
done


mkdir -p vardbs/grch37/strucvar

touch vardbs/grch37/strucvar/{clinvar,dbvar,dgv_gs,dgv,exac,g1k,gnomad_sv,patho_mms}.spec.json
for path in vardbs/grch37/strucvar/{clinvar,dbvar,dgv_gs,dgv,exac,g1k,gnomad_sv,patho_mms}.bed.gz; do
    touch $path
    md5sum $path >$path.md5
done
