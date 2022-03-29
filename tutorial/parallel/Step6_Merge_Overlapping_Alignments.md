[← Back to Step 5](Step5_Alignment_and_Alignment_Trimming.md)


# STEP 6

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step6.png)


## 1) get.consensus.from.alignment.parallel.sh

**Usage**
```
get.consensus.from.alignment.parallel.sh -s <file> -d <directory> -m <positive integer> -b <positive numeric> -t <positive integer> -gnv
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
get.consensus.from.alignment.parallel.sh -s $taxa -d $trimmed -m 1 -b 0.01 -t 20 -gnv
```

## 2) rename.fasta.headers.R

**Usage**
```
rename.fasta.headers.R <file> <string> <BOOLEAN> <BOOLEAN>
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
rename.fasta.headers.R $cons ".all.aln.etr.itr" FALSE FALSE
```

## 3) blast.vs.self.sh

**Usage**
```
blast.vs.self.sh <file>
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
blast.vs.self.sh $cons
```

## 4) find.overlapping.alignments.R

**Usage**
```
find.overlapping.alignments.R <file> <BOOLEAN> <string> <string>
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
find.overlapping.alignments.R $cbase.vs.self.blast.filtered TRUE 'LG_' '_'
```

## 5) align.overlapping.contigs.sh

**Usage**
```
align.overlapping.contigs.sh -l <file> -c <directory> -m <string> -t <positive integer>
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
align.overlapping.contigs.sh -l $overlaps -c $multifasta -m 'localpair' -t 20
```

## 6) filter.merged.alignments.sh

**Usage**
```
filter.merged.alignments.sh -d <directory> -s <numeric fraction>
```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
filter.merged.alignments.sh -d $merged -s $20
```

## 7) replace.overlapping.alignments.R

**Usage**
```
trim.alignment.ends.parallel.sh -s $2 -d $merged -c 0.5 -n 0.25 -t 20 -v
trim.alignments.parallel.sh -s $2 -d $merged -c 0.4 -z 20 -n 0.5 -S 1 -t 20 -v
replace.overlapping.alignments.R <directory> <directory> <file>

```

**Arguments**
```
# Required



# Optional

```

**Depends on**
```

```


**Example**
```
replace.overlapping.alignments.R $trimmed $merged $overlaps
```

## Continue
[➜ Continue to Step 7](Step7_Create_Representative_Reference_Sequences.md)
