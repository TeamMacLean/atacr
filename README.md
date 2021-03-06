[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/TeamMacLean/atacr.svg?branch=master)](https://travis-ci.org/TeamMacLean/atacr)
[![codecov](https://codecov.io/gh/TeamMacLean/atacr/branch/master/graph/badge.svg)](https://codecov.io/gh/TeamMacLean/atacr)
 
-----------------------------------------
 
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.0.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/atacr)](https://cran.r-project.org/package=atacr)
![package_version](https://img.shields.io/badge/Package%20version-0.4.14-orange.svg?style=flat-square)


---------------------------------------
 
![Last-changedate](https://img.shields.io/badge/last%20change-"2018--05--15"-yellowgreen.svg)
<!-- README.md is generated from README.Rmd. Please edit that file -->

atacr
=====

Helps with the analysis of count data from RNA-capture-seq and ATAC-capture-seq experiments. Using BioConductor RangedSummarizedExperiment objects, atacr implements a set of helper functions and quality control plots specific to the analysis of counts of reads in windows across genomes. Especially, atacr is useful for performing sample normalizations and for easily running bootstrap and Bayes factor tests for differentially accessible windows in common reference designs.

Installation
------------

You can install atacr from github with:

``` r
# install.packages("devtools")
devtools::install_github("TeamMacLean/atacr")
```

Documentation
--------------

You can read documentation on the following topics

  1. [Tutorial - A worked example](https://teammaclean.github.io/atacr)
  2. [atacR - General Overview](https://teammaclean.github.io/atacr/atacr.html)
  3. [Loading Data](https://teammaclean.github.io/atacr/loading.html)
  3. [Summaries of Data](https://teammaclean.github.io/atacr/summaries.html)
  4. [Normalising Data](https://teammaclean.github.io/atacr/normalisations.html)
  5. [Differential Windows](https://teammaclean.github.io/atacr/differential_windows.html)
  6. [Subsetting Data](https://teammaclean.github.io/atacr/atacr_which.html)

Quick start:
------------

``` r
library(atacr)
summary(sim_counts)
#> ATAC-seq experiment of 2 treatments in 6 samples
#>  Treatments: control,treatment 
#>  Samples: control_001,control_002,control_003,treatment_001,treatment_002,treatment_003 
#>  Bait regions used: 500 
#>  Total Windows: 1000 
#>  
#>  On/Off target read counts:
#>           sample off_target on_target percent_on_target
#> 1   control_001        312     15160          97.98345
#> 2   control_002        347     14777          97.70563
#> 3   control_003        339     15115          97.80639
#> 4 treatment_001        321     16955          98.14193
#> 5 treatment_002        346     16490          97.94488
#> 6 treatment_003        335     17064          98.07460 
#>  Quantiles: 
#>  $bait_windows
#>     control_001 control_002 control_003 treatment_001 treatment_002
#> 1%        19.99       16.99          19         16.99         16.00
#> 5%        22.00       20.00          22         20.00         19.00
#> 95%       40.00       40.00          39         63.00         65.05
#> 99%       45.00       46.00          44        109.00         89.03
#>     treatment_003
#> 1%          16.00
#> 5%          21.00
#> 95%         61.00
#> 99%        109.06
#> 
#> $non_bait_windows
#>     control_001 control_002 control_003 treatment_001 treatment_002
#> 1%            0           0        0.00             0          0.00
#> 5%            0           0        0.00             0          0.00
#> 95%           3           4        3.05             3          3.05
#> 99%           4           4        4.00             4          4.00
#>     treatment_003
#> 1%              0
#> 5%              0
#> 95%             3
#> 99%             4
#>  
#>  Read depths:
#>           sample off_target on_target
#> 1   control_001      0.624    30.320
#> 2   control_002      0.694    29.554
#> 3   control_003      0.678    30.230
#> 4 treatment_001      0.642    33.910
#> 5 treatment_002      0.692    32.980
#> 6 treatment_003      0.670    34.128
```

``` r
plot(sim_counts)
#> Picking joint bandwidth of 0.0243
#> Picking joint bandwidth of 0.0582
```

![](README-unnamed-chunk-2-1.png)
