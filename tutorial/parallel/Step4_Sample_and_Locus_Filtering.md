[← Back to Step 3](Step3_Orthology_Assessment.md)


# STEP 4

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step4.png)


## 1) [filter.visual.assemblies.sh](https://github.com/scrameri/CaptureAl/wiki/filter.visual.assemblies.sh)

This filters samples and loci (target regions) with poor assembly results, taking taxon groups into account and using filtering thresholds informed by comprehensive visualizations.

**Example**
```
filter.visual.assemblies.sh -s samples.mapfile.txt -t loci_stats.txt -r reference.fasta \
				    -a 0.3 -b 0.5 -c 0.3 -d  1 -e  1 -f 50 \
				    -g 0.2 -h   1 -i   1 \
            -p 0.9
```

## Continue
[➜ Continue to Step 5](Step5_Alignment_and_Alignment_Trimming.md)
