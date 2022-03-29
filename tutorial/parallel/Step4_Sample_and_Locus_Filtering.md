[← Back to Step 3](Step3_Orthology_Assessment.md)


# STEP 4

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step4.png)


## 1) SCRIPT

**Usage**
```
filter.visual.assemblies.sh -s  <file> -t <file> -r <file> \
				    -a <numeric fraction> -b <numeric fraction> -c <numeric fraction> \
            -d <positive integer> -e <positive numeric> -f <positive integer> \
            -g <numeric fraction> -h <positive numeric> -i <positive integer> \
            -p  <numeric fraction>
```

**Arguments**
```
# Required
-s         Path to samples file. Header and tab-separation expected.
          
           Sample IDs must be in the FIRST column. These must match (a subset of) sample names in the mapping 
           stats passed via -t.
          
           Group IDs can be specified in the SECOND column (if not specified, all 
           samples are assumed to constitute one group).
          
           The group ID is used to apply region filtering criteria d-i within all considered groups, to determine regions
           passing the filtering criteria in all groups.
          
           Samples that do not belong to any specified group (second column empty or 'NA') will be displayed in summary
           plots but will not be considerd during region filtering. 
          
           Additional columns are ignored.
          
          
-t         Path to assembly statistics. Header and tab-separation expected.

           Sample IDs must be in the FIRST column. Assembly statistics must be in the following columns as defined in
           filter.visual.coverages.R.
          
           Only assembly statistics of samples passed via -s will be used. A Warning or Stop is issued if there
           are mismatches.


-r         Path to region reference sequences. FASTA format expected. Used to correlate alignment stats with 
           reference sequence lengths and GC content.
          
           Only target regions passed via -t will be considered. A Warning or Stop is issued if there are mismatches.


# Optional
           # The first two filters take absolute thresholds and aim to remove poorly assembled samples or target regions:
-a  [0.3]  minimum fraction of regions with at least one contig in a sample (filters samples)
-b  [0.5]  minimum median contig length relative to target length 
-c  [0.3]  minimum fraction of samples with at least one contig in a region (filters target regions)

           # The next six filters set the thresholds...
-d    [1]  maximum number of non-zero (fragments combined) contigs in a target region
-e    [1]  minimum normalized EXONERATE alignment score
-f   [50]  minimum EXONERATE alignment length
-g  [0.2]  minimum alignment fraction (EXONERATE alignment length divided by target region length)
-h    [1]  minimum raw EXONERATE alignment score
-i    [1]  minimum contig length.

           # ... that need to be met in a specified fraction of samples in each considered taxon group:
-p  [0.9]  minimum fraction of samples in each taxon group that need to pass each filter in order to keep a certain target region.

```

**Depends on**
```
ape
tidyr
```


**Example**
```
filter.visual.assemblies.sh -s samples.mapfile.txt -t loci_stats.txt -r reference.fasta \
				    -a 0.3 -b 0.5 -c 0.3 -d  1 -e  1 -f 50 \
				    -g 0.2 -h   1 -i   1 \
            -p 0.9
```

## Continue
[➜ Continue to Step 5](Step5_Alignment_and_Alignment_Trimming.md)
