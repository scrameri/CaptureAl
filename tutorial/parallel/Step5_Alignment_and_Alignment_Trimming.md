[← Back to Step 4](Step4_Sample_and_Locus_Filtering.md)


# STEP 5

## Overview
![Step.png](https://raw.githubusercontent.com/scrameri/CaptureAl/master/tutorial/CaptureAl_Step5.png)


## 1) create.multifastas.parallel.sh

**Usage**
```
create.multifastas.parallel.sh -s $taxa -l $regions -d $exonerate -t 20
```

**Arguments**
```
# Required



# Optional [DEFAULT]

```

**Depends on**
```

```


**Example**
```

```

## 2) align.multifastas.parallel.sh


**Usage**
```
align.multifastas.parallel.sh -d $multifasta -m 'localpair' -t 20
```

**Arguments**
```
# Required



# Optional [DEFAULT]

```

**Depends on**
```

```


**Example**
```

```

## 3) trim.alignment.ends.parallel.sh


**Usage**
```
trim.alignment.ends.parallel.sh -s $2 -d $mafft -c 0.5 -n 0.25 -t 20 -v
```

**Arguments**
```
# Required
trim.alignment.ends.parallel.sh -s $2 -d $mafft -c 0.5 -n 0.25 -t 20 -v



# Optional [DEFAULT]

```

**Depends on**
```

```


**Example**
```

```


## 4) trim.alignments.parallel.sh


**Usage**
```
trim.alignments.parallel.sh -s $2 -d $endtrimmed -c 0.4 -z 20 -S 1 -n 0.5 -t 20 -v
```

**Arguments**
```
# Required
trim.alignment.ends.parallel.sh -s $2 -d $mafft -c 0.5 -n 0.25 -t 20 -v



# Optional [DEFAULT]

```

**Depends on**
```

```


**Example**
```


## Continue
[➜ Continue to Step 6](Step6_Merge_Overlapping_Alignments.md)
