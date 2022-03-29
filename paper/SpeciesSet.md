## Read quality-trimming and quality-filtering
Raw paired-end reads with the indicated file extensions (`-x` option) located in the `$raw` directory (`-r` option) were quality-trimmed and quality-filtered using trimmomatic version 0.32 (Bolger et al., 2014). Specifically, we used `ILLUMINACLIP` with an adapter sequence file containing NEBNext, TRUSEQ and Illumina adaptor sequences (`-a` option), a seed mismatch of 2, a palindrome clip threshold of 20, a simple clip threshold of 10, a minimum adapter length of 10, while keeping both reads. Leading and trailing bases of each read were removed if the quality was below 5. Sliding window trimming was performed using a window size of 4 and a required average quality of 15. Quality-trimmed reads shorter than 50 bases were removed. The following script executed trimmomatic as specified above, for 20 samples in parallel: 

```
trim.fastq.sh -s samples.dalbergia.12.txt -a illumina.truseq.indexing.adaptors -r $raw -x '_R1.fastq.gz,_R2.fastq.gz' -t 20
```

The quality of raw and trimmed reads was assessed with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc) version 0.11.5.


## Iteration 1 of executed command-line scripts and parameter choices for steps 1–7
The following sections track the executed pipeline scripts and chosen parameters at each analysis step. 

### Step 1: Read mapping
We ran BWA version 0.7.12-r1039 (Li & Durbin, 2009) and BWA-MEM in the `$mapped` directory, using the quality-trimmed and quality-filtered reads with the indicated file extensions (`-e` option, comma-separated string denoting file extensions of forward and reverse reads, respectively) located in the `$trimmed` directory (`-d` option), and the respective reference sequences for each taxon set and iteration (`-r` option). The script outputs reads with a minimum alignment score of `10` (`-T` option), marks secondary hits, and only retains reads with a minimum mapping quality of `10` (`-Q` option) in the final SAM files before compressing them to BAM format. An early filtering of target regions with inadequate coverage across samples prevents time-consuming sequence assembly of target regions that would likely be filtered out in step 4. Computations were performed for all samples specified in `samples.dalbergia.12.txt` (`-s` option) using 4 times 5 threads in parallel (`-t` option) as follows:

```
run.bwamem.sh -s samples.dalbergia.12.txt -r Cajanus_cajan_6555reg.fasta -e .trim1.fastq.gz,.trim2.fastq.gz -T 10 -Q 10 -d $trimmed -t 4
```

We performed coverage analysis on the BAM files filtered for mapping quality equal or above 10 (`-Q` option) for each sample, and wrote all coverage results to one file as follows:

```
get.coverage.stats.sh -s samples.dalbergia.12.txt -Q 10 -t 20

collect.coverage.stats.R samples.dalbergia.12.txt 10
```

This produced the file **`coverage_stats.txt`**, which was used to perform target region filtering. We implemented seven filtering criteria to identify target regions with adequate average coverage across the taxon groups specified in `mapfile.dalbergia.12.txt`.

The first two filters take *absolute* thresholds and aim to remove poorly sequenced samples or target regions:
- `-a` = `min.pregion` = minimum fraction of regions with at least one mapped read in a sample (filters samples)
- `-b` = `min.ptaxa`   = minimum fraction of samples with at least one mapped read in a region (filters target regions)

The next four filters take thresholds...
- `-c` = `min.len`     = minimum BWA-MEM alignment length
- `-d` = `min.cov`     = minimum average coverage in the aligned region
- `-e` = `max.cov`     = maximum average coverage in the aligned region
- `-f` = `min.ratio`   = minimum alignment fraction (BWA-MEM alignment length divided by target region length)

...that need to be met in a specified *fraction* of samples in each considered taxon group:
- `-p` = `min.frac`    = minimum fraction of samples in each taxon group that need to pass each filter in order to keep a certain target region

```
filter.visual.coverages.sh -s mapfile.dalbergia.12.txt -t coverage_stats.txt -r Cajanus_cajan_6555reg.fasta -a 0.2 -b 0.4 -c 1 -d 8 -e 1000 -f 0 -p 0.4
```

This script visualized the coverage statistics as violin plots and heatmaps (see Figure S6 for results of the second iteration), and saved a list of kept samples (`coverage_stats-taxa-*.txt`) as well as a list of kept target regions (`coverage_stats-regions-*.txt`) for sequence assembly.


### Step 2: Sequence assembly
We extracted read pairs from quality-filtered and quality-trimmed reads located in the `$trimmed` directory (`-d` option). The `-s` and `-l` parameters are used to pass the list of samples and loci (target regions) to be processed in parallel, respectively. This step was carried out on a local scratch (`$scratch` directory) using 20 parallel threads (`-t` option) and produced the `$extractedreads` directory. At least one of the two reads per extracted read pair mapped to a retained target region with a minimum mapping quality of 10 (`-Q` option):  

```
s=coverage_stats-taxa-0.2.txt
l=coverage_stats-regions-0.2-0.4-1-8-1000-0-0.4.txt

cd $scratch

extract.readpairs.sh -s $s -l $l -d $trimmed -m $mapped -Q 10 -t 20
```

We assembled the extracted reads located in the `$extractedreads` directory (`-r` option) into consensus contigs (contigs hereafter) separately for each sample and retained region using DIPSPADES (SPADES version 3.6.0) in `assembly-only` and `careful` mode, with an automatic coverage cutoff. This step was carried out on a local scratch (`$assemblies` directory) using 20 parallel threads (`-t` option):

```
run.dipspades.sh -s $s -r $extractedreads -t 20
```


### Step 3: Orthology assessment
We ran EXONERATE version 2.2 for each sample and each retained target region (`-l` option) with the `affine:local` and `exhaustive` options, using the contigs located in the `$assemblies` directory (`-d` option) as query sequences and the target regions (`-r` option) as target sequences. We stored alignment statistics of all consensus contigs that aligned to the same target region in the `exonerate` directory, but limited the report to the best alignment per contig as follows:

```
select.best.contigs.per.locus.sh -s $s -l $l -r Cajanus_cajan_6555reg.fasta -d $assemblies -t 20
```

Contigs with a target alignment length of at least the specified threshold (`-a` option) and a normalized alignment score (defined as the raw EXONERATE alignment score divided by the target alignment length) of at least the specified threshold (`-c` option) were considered as potentially homologous and retained. If more than one contig met these requirements, and if none of these contigs physically overlapped based on the alignment statistics, the best-matching contig was combined with the additional contig(s) using an appropriate spacer and by taking the directionality into account as follows:

```
combine.contigs.parallel.sh -s $s -d $exonerate -a 1 -c 2 -t 20
```

We collected the EXONERATE statistics of each sample and plotted the number of contigs per target region for the different taxon groups as follows:

```
collect.exonerate.stats.R $s $exonerate
```

This produced the file `loci_contignumbers.txt`, which was visualized using the following R script, as well as the file **`loci_stats.txt`**, which contained the EXONERATE statistics used to perform sample and target region filtering in step 4.

```
plot.contig.numbers.R loci_contignumbers.txt mapfile.dalbergia.12.txt
```

This produced Figure S8c (results of the second iteration shown).


### Step 4: sample and region filtering
We implemented nine filtering criteria to identify target regions with adequate assembly quality across the taxon groups specified in `mapfile.dalbergia.12.txt`.

The first two filters take absolute thresholds and aim to remove poorly assembled samples or target regions:
- `-a` = `min.pregion`        = minimum fraction of regions with at least one contig in a sample (filters samples)
- `-b` = `min.pcontig`        = minimum median contig length relative to target length 
- `-c` = `min.ptaxa`          = minimum fraction of samples with at least one contig in a region (filters target regions)

The next six filters set the thresholds...
- `-d` = `max.ncontigs`       = maximum number of non-zero (fragments combined) contigs in a target region
- `-e` = `min.bestscore.norm` = minimum normalized EXONERATE alignment score
- `-f` = `min.taln`           = minimum EXONERATE alignment length
- `-g` = `min.tfrac`          = minimum alignment fraction (EXONERATE alignment length divided by target region length)
- `-h` = `min.bestscore`      = minimum raw EXONERATE alignment score
- `-i` = `min.bestlength`     = minimum contig length.

... that need to be met in a specified fraction of samples in each considered taxon group:
- `-p` = `min.frac`           = minimum fraction of samples in each taxon group that need to pass each filter in order to keep a certain target region.

```
filter.visual.assemblies.sh -s mapfile.dalbergia.12.txt -t loci_stats.txt -r Cajanus_cajan_6555reg.fasta -a 0.2 -b 0 -c 0.5 -d 2 -e 2 -f 80 -g 0 -h 1 -i 1 -p 0.5
```

This script visualized the assembly statistics as violin plots and heatmaps (see panels a, b and d of Figure S8 for results of the second iteration), and saved a list of kept samples (`taxa=taxa_kept-*.txt`) as well as a list of kept target regions (`regions_kept-*.txt`) for sequence alignment.


### Step 5: Target region alignment and alignment trimming
We generated multifasta files in the `$multifasta directory for all retained target regions, containing all retained contigs and samples as follows:

```
taxa=taxa_kept-0.2.txt
regions=regions_kept-0.2-0-0.5-2-2-80-0.5.txt

create.multifastas.parallel.sh -s $taxa -l $regions -d $exonerate -t 20
```

We generated alignments in the $mafft directory using MAFFT version 7.123b, the ‘localpair’ and ‘adjustdirection’ flags, and 1000 maximum iterations:

```
align.multifastas.parallel.sh -d $multifasta -m 'localpair' -t 20
```

Raw alignments were trimmed at both ends until an alignment site had nucleotides in at least 50% of aligned sequences (`-c` option) and a maximum nucleotide diversity (i.e., the sum of the number of base differences between sequence pairs divided by the number of comparisons) of 0.25 (`-n` option). The `-v` flag triggered visualization of the alignment end trimming procedure. Trimmed alignments were written to the `$endtrimmed` directory as follows:

```
trim.alignment.ends.parallel.sh -s $taxa -d $mafft -c 0.5 -n 0.25 -t 20 -v
```

Internal trimming was carried out by first removing any alignment site with nucleotides in less than 40% of aligned sequences (`-c` option). Potential mis-assemblies or mis-alignments in each sequence were resolved using a sliding window approach with window size 20 (`-z` option) and step size 1 (`-S` option). Specifically, we trimmed windows at contig ends if more than 50% of the nucleotides in the conserved part of the window deviated from the alignment consensus (`-n` option). The script defines a conserved part of each window as the alignment sites with nucleotides in at least 20% of samples, and where the frequencies of minor alleles are all below 30% without considering gaps. After window-based trimming, the script also removes sites with sequence data for less than the specified fraction of aligned sequences (`-c` option) again. The `-v` flag triggered visualization of the internal trimming procedure with sliding window approach. Trimmed alignments were written to the `$trimmed` directory as follows:

```
trim.alignments.parallel.sh -s $taxa -d $endtrimmed -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -v
```

### Step 6: Merge overlapping alignments
We calculated a consensus sequence for each end-trimmed and internally trimmed alignment located in the directory `$trimmed`, using a minimum allele frequency of 1 (`-m` option) to call IUPAC ambiguity and a minimum base frequency of 0.01 (`-b` option) to return a consensus instead of a gap. This parameter combination ensured that the most frequent allele was called at each alignment site rather than IUPAC ambiguity codes or gaps. If two alleles were equally frequent at any alignment site, one was randomly sampled to represent the consensus. The `-g` flag ensured that gaps were removed from the final consensus sequence, and the `-n` flag ensured that completely ambiguous consensus bases (Ns) were removed from the final consensus sequence. The `-v` flag triggered visualization of the consensus calculation:

```
get.consensus.from.alignment.parallel.sh -s $taxa -d $trimmed -m 1 -b 0.01 -t 20 -gnv
```

We renamed the sequence names of alignment consensus sequences stored in `$cons` to dispose of the suffix added during alignment and trimming before identifying the best non-reciprocal BLAST+ hits between alignment consensus sequences as follows:

```
cons=${trimmed}-cons.fasta
rename.fasta.headers.R $cons ".all.aln.etr.itr" FALSE FALSE
blast.vs.self.sh $cons
```

We then identified BLAST+ hits at alignment ends and stored a list of physically overlapping alignments (names of overlapping target regions on the same line) as follows:

```
cbase=$(basename $cons .fasta)
find.overlapping.alignments.R $cbase.vs.self.blast.filtered TRUE 'LG_' '_'
```

Arguments 2–4 limit the identification of overlapping alignments to target regions of the same linkage group, which is identified via user-specified strings surrounding a linkage group identifier present in target region names. We then aligned all contigs of up to five physically overlapping target regions using the same alignment algorithm as before. Merged alignments were written to the $merged directory as follows (visualization was triggered by default):

```
overlaps=$cbase.list
align.overlapping.contigs.sh -l $overlaps -c $multifasta -m 'localpair' -t 20
```

In cases where contigs of the same sample overlapped with a mismatch, only the base with higher frequency at that alignment site was considered. A success score of each merging procedure was computed based on the number of mismatches in overlapping contigs of the same individuals relative to the total number of bases in the alignment. The success score amounted to 1 if there were no mismatches in any individual. We discarded any merged alignment with a score smaller than 0.85 (subfamily set) or 0.9 (species set), passed to the `-s` option.

```
filter.merged.alignments.sh -d $merged -s 0.85
```

Unsuccessfully merged alignment sets were visually inspected to identify whether some subsets of alignments sufficiently overlapped to allow for merging. The manually selected alignments were merged again and combined with the automatically merged alignments if they showed a sufficient success score. All successfully merged alignments were then trimmed as before and used as replacements for overlapping alignments.

```
trim.alignment.ends.parallel.sh -s $taxa -d $merged -c 0.5 -n 0.25 -t 20 -v
trim.alignments.parallel.sh -s $taxa -d $merged -c 0.4 -z 20 -n 0.5 -S 1 -t 20 -v
replace.overlapping.alignments.R $trimmed $merged $overlaps
```

### Step 7: Create representative reference sequences
Sets of reference consensus sequences for different taxon groups were generated, combined, aligned, and a group consensus was derived as follows:

```
get.group.consensus.sh -s mapfile.dalbergia.12.txt -d $trimmed -m 1 -b 0.01 -z ".all.aln.etr.itr.cons" -t 20 -gnv
```

The resulting FASTA file $newref was renamed according to the taxon set and number of remaining target regions:

```
rename.fasta.headers.R $newref ".cons.aln" FALSE FALSE
mv $newref Dalbergia_iter1_2468reg.fasta
```


## Iteration 2 of executed command-line scripts and parameter choices for steps 1–7
Iteration 1 was applied to a subset of 12 samples, which were representative of all defined taxon groups. We applied a second iteration of read mapping, assembly, and alignment to all 110 samples of the subfamily set, using more stringent filtering parameters as follows (see Iteration 1 for explanations):


### Step 1: Read mapping
```
run.bwamem.sh -s samples.dalbergia.txt -r Dalbergia_iter1_3736reg.fasta -e .trim1.fastq.gz,.trim2.fastq.gz -T 10 -Q 10 -d $trimmed -t 4
```
```
get.coverage.stats.sh -s samples.dalbergia.txt -Q 10 -t 20

collect.coverage.stats.R samples.dalbergia.txt 10
```
```
filter.visual.coverages.sh -s mapfile.dalbergia.txt -t coverage_stats.txt -r Dalbergia_iter1_3736reg.fasta -a 0.2 -b 0.7 -c 1 -d 8 -e 1000 -f 0 -p 0.7
```

### Step 2: Sequence assembly
```
s=coverage_stats-taxa-0.2.txt
l=coverage_stats-regions-0.2-0.7-1-8-1000-0-0.7.txt

cd $scratch

extract.readpairs.sh -s $s -l $l -d $trimmed -m $mapped -Q 10 -t 20
```
```
run.dipspades.sh -s $s -r $extractedreads -t 20
```


### Step 3: Orthology assessment
```
select.best.contigs.per.locus.sh -s $s -l $l -r Dalbergia_iter1_3736reg.fasta -d $assemblies -t 20
```
```
combine.contigs.parallel.sh -s $s -d $exonerate -a 1 -c 2 -t 20
```
```
collect.exonerate.stats.R $s $exonerate
```
```
plot.contig.numbers.R loci_contignumbers.txt mapfile.dalbergia.txt
```


### Step 4: sample and region filtering
```
filter.visual.assemblies.sh -s mapfile.dalbergia.txt -t loci_stats.txt -r Dalbergia_iter1_3736reg.fasta -a 0.2 -b 0 -c 0.75 -d 2 -e 2 -f 80 -g 0 -h 1 -i 1 -p 0.5
```


### Step 5: Target region alignment and alignment trimming
```
taxa=taxa_kept-0.2.txt
regions=regions_kept-0.2-0-0.75-2-2-80-0.5.txt

create.multifastas.parallel.sh -s $taxa -l $regions -d $exonerate -t 20
```
```
align.multifastas.parallel.sh -d $multifasta -m 'localpair' -t 20
```
```
trim.alignment.ends.parallel.sh -s $taxa -d $mafft -c 0.5 -n 0.25 -t 20 -v
```
```
trim.alignments.parallel.sh -s $taxa -d $endtrimmed -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -v
```


### Step 6: Merge overlapping alignments
```
get.consensus.from.alignment.parallel.sh -s $taxa -d $trimmed -m 1 -b 0.01 -t 20 -gnv
```
```
cons=${trimmed}-cons.fasta
rename.fasta.headers.R $cons ".all.aln.etr.itr" FALSE FALSE
blast.vs.self.sh $cons
```
```
cbase=$(basename $cons .fasta)
find.overlapping.alignments.R $cbase.vs.self.blast.filtered TRUE 'LG_' '_'
```
```
overlaps=$cbase.list
align.overlapping.contigs.sh -l $overlaps -c $multifasta -m 'localpair' -t 20
```
```
filter.merged.alignments.sh -d $merged -s 0.85
```
```
trim.alignment.ends.parallel.sh -s $taxa -d $merged -c 0.5 -n 0.25 -t 20 -v
trim.alignments.parallel.sh -s $taxa -d $merged -c 0.4 -z 20 -n 0.5 -S 1 -t 20 -v
replace.overlapping.alignments.R $trimmed $merged $overlaps
```


### Step 7: Create representative reference sequences
```
get.group.consensus.sh -s mapfile.dalbergia.txt -d $trimmed -m 1 -b 0.01 -z ".all.aln.etr.itr.cons" -t 20 -gnv
```
```
rename.fasta.headers.R $newref ".cons.aln" FALSE FALSE
mv $newref consDalbergia_4c_2396.fasta 
```
