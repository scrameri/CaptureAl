# Handling Target Enrichment Probe Sequences

## Merge physically neighboring reference sequences based on reference genome
- [merge_refseqs.R](merge_refseqs.R)
  - Combine two phylocally close sequences into a single contiguous sequence using a reference genome to fill the gap.

## Filter BLAST+ hits for reciprocal best hits with minimum alignment length
- [filter_BLAST_Angiosperms353.R](filter_BLAST_Angiosperms353.R)
  - Filter raw BLAST+ output for reciprocal best hits with minimum alignment length.
