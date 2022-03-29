# Supplementary Methods
[CaptureAl](https://github.com/scrameri/CaptureAl) accompanies a manuscript in review (Crameri et al. submitted), available [here](https://www.biorxiv.org/content/10.1101/2021.12.07.471551v1). This repository contains some further documentation on Supplementary Methods.

## Merge physically neighboring reference sequences
- [merge_refseqs.R](ProbeSets/merge_refseqs.R)
  - Combine two phylocally close sequences into a single contiguous sequence using a reference genome to fill the gap.

## Complete workflow of executed CaptureAl pipeline steps
- [Subfamily set](SubfamilySet.md)
- [Species set](SpeciesSet.md)

## Create new probe sets
- [extract_probes.R](ProbeSets/extract_probes.R)
  - Extract tiled target enrichment probes from a list of target regions

## Filter BLAST+ hits for reciprocal best hits with minimum alignment length
- [filter_BLAST_Angiosperms353.R](ProbeSets/filter_BLAST_Angiosperms353.R)
  - Filter raw BLAST+ output for reciprocal best hits with minimum alignment length.
