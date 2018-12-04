Gene Expression Analysis (of SNPs)
================
Raghav Chanchani
11/15/2018

Asthma in People of Mexican Ancestry in Los Angeles, CA
-------------------------------------------------------

``` r
asthma <- read.csv("asthma.csv")
```

Focusing on the second column which shows the subjects' genotype... The genotype is formatted as the diploid alleles separated by "|". % different alleles are shown below.

``` r
genotypes <- round(table(asthma[,2]) / nrow(asthma) * 100, 2)
```

There are 34.38% people with the A|A genotype.

Interpreting Base Qualities
---------------------------

``` r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Expression Data Analysis
------------------------

``` r
table(expr$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
inds.ag <- expr$geno =="AG"
summary(expr$exp[inds.ag])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 

``` r
inds.gg <- expr$geno =="GG"
summary(expr$exp[inds.gg])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ## 

``` r
boxplot(exp ~ geno, data=expr) 
```

![](class14_files/figure-markdown_github/unnamed-chunk-10-1.png)
