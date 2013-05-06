Whole Genome Bisulfite Sequencing
=================================

Overview
--------
This example demonstrates the use of our smoothing approach for improving the
signal-to-noise ratio for Whole Genome Bisulfite Sequencing (WGBS) methylation
data.

(...methylation background...)

WGBS data is inherently noisy, and as such the signatures of interesting
features such as CpG islands is often obscured by numerous low coverage
reads.

The data used for this analysis comes from a study by Hansen et al. (2011)^1,
and consists of three components:

1. **cpgsites** - Genomic location (offset) of CpG sites in the dataset.
2. **methylation** - Number of methylation reads at a given CpG site.
3. **coverage** - Total number of reads at a given site.

The first variable (cpgsites) is a mx1 integer vector, while the later two
variables are each mxn matrices, each with $n = 6$ separate biological samples.

The smoothing method smooths the methylation values across time, weighting
the contribution of each site by the ratio of the number of methylated reads
to total reads at the site.

Analysis
--------
First let's load sample data and take a look at what we have:

```r
source("../smooth.r")
load(file = "../data/sample.rda")

# Methylation reads at each site
head(methylation)
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    1    5    4    9    1    2
## [2,]    3    6    4    6    1    3
## [3,]    3    3    3    2    2    1
## [4,]    3    2    3    4    1    3
## [5,]    3    2    2    2    1    2
## [6,]    2    1    1    2    4    0
```

```r
summary(methylation)
```

```
##        V1              V2              V3              V4       
##  Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
##  1st Qu.:    1   1st Qu.:    1   1st Qu.:    2   1st Qu.:    1  
##  Median :    4   Median :    3   Median :    4   Median :    3  
##  Mean   :    4   Mean   :    3   Mean   :    4   Mean   :    4  
##  3rd Qu.:    6   3rd Qu.:    5   3rd Qu.:    6   3rd Qu.:    6  
##  Max.   :16306   Max.   :15501   Max.   :16506   Max.   :13777  
##        V5              V6       
##  Min.   :    0   Min.   :    0  
##  1st Qu.:    1   1st Qu.:    1  
##  Median :    3   Median :    2  
##  Mean   :    4   Mean   :    3  
##  3rd Qu.:    6   3rd Qu.:    5  
##  Max.   :13826   Max.   :14199
```

```r

# Total number of reads at each site
head(coverage)
```

```
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    4    5    4    9    1    2
## [2,]    4    6    4    7    1    3
## [3,]    4    3    3    2    2    1
## [4,]    4    2    4    4    1    3
## [5,]    3    2    2    3    2    3
## [6,]    2    2    2    5    4    0
```

```r
summary(coverage)
```

```
##        V1              V2              V3              V4       
##  Min.   :    0   Min.   :    0   Min.   :    0   Min.   :    0  
##  1st Qu.:    3   1st Qu.:    3   1st Qu.:    3   1st Qu.:    3  
##  Median :    5   Median :    5   Median :    6   Median :    6  
##  Mean   :    6   Mean   :    5   Mean   :    6   Mean   :    6  
##  3rd Qu.:    8   3rd Qu.:    7   3rd Qu.:    8   3rd Qu.:    8  
##  Max.   :17793   Max.   :16594   Max.   :18090   Max.   :15064  
##        V5              V6       
##  Min.   :    0   Min.   :    0  
##  1st Qu.:    3   1st Qu.:    3  
##  Median :    5   Median :    5  
##  Mean   :    6   Mean   :    5  
##  3rd Qu.:    8   3rd Qu.:    7  
##  Max.   :15222   Max.   :16981
```

```r

# Methylation locations (offset in Chromosome)
head(cpgsites)
```

```
##       [,1]
## [1,] 10469
## [2,] 10471
## [3,] 10484
## [4,] 10489
## [5,] 10493
## [6,] 13079
```

```r
range(cpgsites)
```

```
## [1]     10469 249239887
```


Let's now apply our smoothing method and compare the results.

```r
result = fast.smooth(cpgsites, methylation, coverage)

# down-sampled data
head(result$small)
```

```
##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
## [1,] 0.8392 0.9307 0.9100 0.8713 0.9151 0.9475
## [2,] 0.8169 0.8850 0.8532 0.7859 0.9072 0.8614
## [3,] 0.7191 0.5577 0.6368 0.3119 0.8736 0.4913
## [4,] 0.6491 0.5264 0.6013 0.3846 0.8520 0.5455
## [5,] 0.2500 0.0000 0.7000 0.3333 0.4286 0.5000
## [6,] 1.0000 0.0000 1.0000 0.3333 0.5000 0.5000
```

```r

# after interpolation back up to it's full size
head(result$y)
```

```
##        [,1]   [,2]   [,3]   [,4]   [,5]   [,6]
## [1,] 0.8299 0.9115 0.8862 0.8354 0.9118 0.9113
## [2,] 0.8299 0.9115 0.8862 0.8354 0.9118 0.9113
## [3,] 0.8299 0.9115 0.8861 0.8354 0.9118 0.9113
## [4,] 0.8299 0.9115 0.8861 0.8354 0.9118 0.9113
## [5,] 0.8298 0.9115 0.8861 0.8353 0.9118 0.9113
## [6,] 0.8275 0.9067 0.8802 0.8265 0.9109 0.9023
```


Finally, let's create some simple plots of our results.


```r
library(ggplot2)

# plot down-sampled data
dat = stack(as.data.frame(result$small))
scale_factor = (max(cpgsites) - min(cpgsites))/nrow(result$small)
dat$x <- rep(seq_len(nrow(result$small)), ncol(result$small)) * scale_factor
ggplot(dat, aes(x, values)) + geom_line(aes(colour = ind))
```

![plot of chunk visualization](figure/visualization1.png) 

```r

# full-size
dat = stack(as.data.frame(result$y))
dat$x <- rep(seq_len(nrow(result$y)), ncol(result$y))
ggplot(dat, aes(x, values)) + geom_line(aes(colour = ind))
```

```
## Warning: Removed 126 rows containing missing values (geom_path).
```

![plot of chunk visualization](figure/visualization2.png) 


Done!

References
----------
- Kasper Daniel Hansen, Winston Timp, Héctor Corrada Bravo, Sarven Sabunciyan, Benjamin Langmead, Oliver G McDonald, Bo Wen, Hao Wu, Yun Liu, Dinh Diep, Eirikur Briem, Kun Zhang, Rafael A Irizarry, Andrew P Feinberg,   (2011) Increased Methylation Variation in Epigenetic Domains Across Cancer Types.  *Nature Genetics*  **43**  [10.1038/ng.865](http://dx.doi.org/10.1038/ng.865)

