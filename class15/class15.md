Genome Informatics II
================
Raghav Chanchani
11/20/2018

RNA-Seq Analysis
----------------

First step is to read our metadata (countData and colData) files that we will use with DESeq.

``` r
counts <- read.csv("data/airway_scaledcounts.csv", stringsAsFactors = FALSE)
metadata <- read.csv("data/airway_metadata.csv", stringsAsFactors = FALSE)
```

Examine the metadata to find **control** and **treated** columns (cell-lines).

``` r
# finding data for control
control <- metadata[metadata$dex == "control", ]
control
```

    ##           id     dex celltype     geo_id
    ## 1 SRR1039508 control   N61311 GSM1275862
    ## 3 SRR1039512 control  N052611 GSM1275866
    ## 5 SRR1039516 control  N080611 GSM1275870
    ## 7 SRR1039520 control  N061011 GSM1275874

``` r
# average of the counts for control
control.mean <- rowSums( counts[ ,control$id])/nrow(control)
# adding back gene names
names(control.mean) <- counts$ensgene
```

Doing the same thing for the treated trials...

``` r
# finding data for treated
treated <- metadata[metadata$dex == "treated", ]
treated
```

    ##           id     dex celltype     geo_id
    ## 2 SRR1039509 treated   N61311 GSM1275863
    ## 4 SRR1039513 treated  N052611 GSM1275867
    ## 6 SRR1039517 treated  N080611 GSM1275871
    ## 8 SRR1039521 treated  N061011 GSM1275875

``` r
# average of counts for treated
treated.mean <- rowSums( counts[ ,treated$id])
# adding back gene names
names(treated.mean) <- counts$ensgene
```
