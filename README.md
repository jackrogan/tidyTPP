
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyTPP

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/jackrogan/tidyTPP/graph/badge.svg?token=OFYPL51CV2)](https://app.codecov.io/gh/jackrogan/tidyTPP)
[![R-CMD-check](https://github.com/jackrogan/tidyTPP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jackrogan/tidyTPP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Introduction

tidyTPP provides an analysis pipeline for *Thermal Protein Profiling
(TPP)* proteomics data, with functions for reading in data from protein
quantification software, normalising data, analysis, and hit-finding.

The goal of tidyTPP is to provide a range of functions, each of which
requires data in a straightforward data frame or
[tibble](https://tibble.tidyverse.org/) and returns long-format data in
a [tibble](https://tibble.tidyverse.org/) that are readily understood
and easily used for further analysis. The style of data is inspired by
the [tidyverse](https://www.tidyverse.org/), and tidyTPP makes use of
several tidyverse packages - in particular
[tibble](https://tibble.tidyverse.org/) and
[ggplot2](https://ggplot2.tidyverse.org/).

### Currently supports:

Importing:

- Importing from Thermo Proteome Discoverer

- Importing from Spectronaut

- Importing with custom settings

Normalisation and Analysis:

- Single- and multi-threaded fitting of temperature-dependent protein
  melting curves *(TPP-TR)*

- Significance scores by melting point difference
  ($\Delta T_m$)<sup>[1](#ref-Savitski_2014)</sup>

- Significance scores by *non-pararametric analysis of response curves
  (NPARC)*<sup>[2](#ref-Childs_2019)</sup>

Hit-finding

- Identifying hits by $\Delta T_m$ p-value

- Identifying hits by *NPARC* F-score and p-value

Plotting

- Plotting *TPP-TR* melting-curves with
  [ggplot2](https://ggplot2.tidyverse.org/)

## Installation

You can install the development version of tidyTPP from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jackrogan/tidyTPP")
```

## Example: Pre-built pipeline

A version of the entire pipeline, with defaults to match common uses, is
available with *apply_TPP_pipeline*. This executes the following
sequence:

``` r
import_TPP() |>
  normalise_TPP() |>
  analyse_TPP() |>
  export_TPP() |>
  get_TPP_hits()
```

The function uses the following defaults in addition to the default
behaviour of each included function:

- Import supported proteomics quantification data formats

- Calculate *NPARC* and $\Delta T_m$ significance scores

- Report all hits with:

  - Adjusted *NPARC*-derived p-value \< 0.05

  - $\Delta T_m$ in the same direction

  - $\Delta T_m$ *vs.* control \> $\Delta T_m$ between controls

  - Full results table and identified hits are exported as *.xlsx* files
    with input file names modified with “\_Results.xlsx” and
    “\_Hits.xlsx” respectively, in the same location as the initial file
    input.

Usage:

``` r
# Example pipeline call
apply_TPP_pipeline(
  datafile = "path_to_data_input.csv",
  config = "path_to_config_file.csv",
  path = "[optional]_path_to_shared_dir",
  import_format = "format_name",
  to_plot = FALSE,
  max_cores = 4)
```

The character argument *import_format* defines which of the inbuilt
import functions to use to read in data files. Currently:

- “Spectronaut”, “SN” - imports from peptide tables exported from
  *Spectronaut* DIA quantification reports.

- “ProteomeDiscoverer”, “PD” - imports from protein tables exported from
  *Thermo Proteome Discoverer* DDA quantification results.

- Default “Spectronaut”.

The (boolean) argument *to_plot* defines whether to show automated plots
across all methods; if *TRUE*, normalisation, NPARC score distribution,
and hit melting-point curves will all be plotted.

- Default *FALSE*.

The (integer) argument *max_cores* is passed to curve-fitting methods
and defines maximum parallel cores to use for these operations.

- Default *4*.

``` r
# Pipeline example: 20-protein test data
library(tidyTPP)

twenty_prot_report <-
    system.file("extdata", "20_protein_peptide_report.csv", package = "tidyTPP")
experiment_config <-
    system.file("extdata", "20_protein_config.csv", package = "tidyTPP")

# Apply pipeline
TPP_hits <- 
  apply_TPP_pipeline(datafile = twenty_prot_report,
                    config = experiment_config,
                    import_format = "spectronaut",
                    to_plot = TRUE,
                    max_cores = 4)
```

Result:

``` r
TPP_hits
#> # A tibble: 2 × 10
#>   Protein_ID Condition Comparison         F_scaled p_adj_NPARC max_adj_pvalue
#>   <chr>      <chr>     <chr>                 <dbl>       <dbl>          <dbl>
#> 1 Protein_I  Treated   Treated_vs_Control     23.2   0.0000901      0.0000141
#> 2 Protein_C  Treated   Treated_vs_Control     14.1   0.00118        0.0981   
#> # ℹ 4 more variables: mean_melt_point <dbl>, mean_control_melt_point <dbl>,
#> #   mean_diff_melt_point <dbl>, mean_control_diff_melt_point <dbl>
```

<img src="man/figures/README-apply_pipeline-4.png"
style="width:100.0%" />

## Example: Function detail

This example will walk through the main functions using 20-protein
example data

### Import

*import\_* functions read two files - quantification results output in a
program-specific format, and a configuration file, defining conditions,
replicates and temperatures, in the format shown.

``` r
# 20-protein example data
twenty_prot_report <-
    system.file("extdata", "20_protein_peptide_report.csv", package = "tidyTPP")
experiment_config <-
    system.file("extdata", "20_protein_config.csv", package = "tidyTPP")
experiment_dir <- system.file("extdata", package = "tidyTPP")
  
# Config file contents
read.csv(experiment_config)[1:10,]
#>    Experiment Condition Replicate Temp
#> 1           1   Control         1   37
#> 2           2   Control         1   41
#> 3           3   Control         1   44
#> 4           4   Control         1   47
#> 5           5   Control         1   50
#> 6           6   Control         1   53
#> 7           7   Control         1   56
#> 8           8   Control         1   59
#> 9           9   Control         1   63
#> 10         10   Control         1   67

# Import using spectronaut import format
twenty_prot_quan_data <- 
  import_spectronaut(datafile = "20_protein_peptide_report.csv",
                     config = "20_protein_config.csv",
                     path = experiment_dir)
#> 
#> --------------------
#> TPP Data Import
#> --------------------
#> Read in:
#> 20_protein_peptide_report.csv 
#> 20_protein_config.csv 
#> --------------------
#> Pivoting to long table...
#> Transforming experiment names...
#> Matching to experiment config data...
#> Finding relative quantity values...
#> Data imported.
#> Found 20 proteins.
#> --------------------

# Resulting tibble
twenty_prot_quan_data
#> # A tibble: 800 × 8
#>    Protein_ID Pep_N Match_N Condition Replicate  Temp rel_quantity raw_quantity
#>    <chr>      <dbl>   <dbl> <chr>     <chr>     <int>        <dbl>        <dbl>
#>  1 Protein_A     29      42 Control   01           37       1          1401047.
#>  2 Protein_A     29      42 Control   01           41       1.05       1464157.
#>  3 Protein_A     29      42 Control   01           44       1.06       1485400.
#>  4 Protein_A     29      42 Control   01           47       1.12       1568061.
#>  5 Protein_A     29      42 Control   01           50       1.35       1889459 
#>  6 Protein_A     29      42 Control   01           53       1.29       1801219.
#>  7 Protein_A     29      42 Control   01           56       0.450       630379.
#>  8 Protein_A     29      42 Control   01           59       0.116       162054.
#>  9 Protein_A     29      42 Control   01           63       0.0945      132437.
#> 10 Protein_A     29      42 Control   01           67       0.166       232101.
#> # ℹ 790 more rows
```

Data at this stage can be plotted to directly observe relative but *not
yet normalised* protein curves.

``` r
plot_melt(twenty_prot_quan_data)
```

<img src="man/figures/README-example_import_plot-1.png" width="100%" />

### Normalisation

*normalise_TPP* transforms relative TPP-TR relative intensity data,
normalising against fitted median melting curves, as described by
Savitsky *et al.* 2014.<sup>[1](#ref-Savitski_2014)</sup>

``` r
# Normalise twenty-protein data, with visualisation
twenty_prot_normalised <- 
  normalise_TPP(TPP_tbl = twenty_prot_quan_data,
                to_plot = TRUE)
#> --------------------
#> TPP Normalisation
#> --------------------
#> Quality Criteria:
#>     col lower upper
#> 1 Pep_N     2   Inf
#> 
#> jointP contains 20 Proteins.
#> 
#> normP criteria:
#>   Temp lower upper
#> 1   56   0.4   0.6
#> 2   63  -Inf   0.3
#> 3   67  -Inf   0.2
#> 
#> normP contains 6 Proteins.
#> --------------------
#> Fit melting curve to normP medians:
#> 
#> Fitting 4 protein curves...
#> Estimated total process time: 2.24 s
#> |==        | 1 of 4|=====     | 2 of 4|=======   | 3 of 4|==========| 4 of 4
#> 4 of 4 fitted successfully.
#> 
#> Total elapsed time: 0.03 s
#> 
#> Best fitted normP median curve:
#>   Condition Replicate  R_sq
#> 1   Control        01 0.995
#> --------------------
#> 20 proteins normalised
#> --------------------
```

<img src="man/figures/README-example_normalisation-1.png"
style="width:100.0%" />

Resulting in the normalised data:

``` r
# Plotted melting data points
plot_melt(twenty_prot_normalised)
#> # A tibble: 800 × 9
#>    Condition Replicate  Temp Protein_ID Pep_N Match_N rel_quantity raw_quantity
#>    <chr>     <chr>     <int> <chr>      <dbl>   <dbl>        <dbl>        <dbl>
#>  1 Control   01           37 Protein_A     29      42            1     1401047.
#>  2 Control   01           37 Protein_B     23      39            1    18637136 
#>  3 Control   01           37 Protein_C     36      62            1     7913466.
#>  4 Control   01           37 Protein_D     76     108            1     1695097.
#>  5 Control   01           37 Protein_T     19      29            1     5209475 
#>  6 Control   01           37 Protein_Q     14      17            1     2789951.
#>  7 Control   01           37 Protein_N     19      35            1      997140.
#>  8 Control   01           37 Protein_K     13      15            1     6945545 
#>  9 Control   01           37 Protein_H     14      26            1     2126099 
#> 10 Control   01           37 Protein_E     10      12            1      351740.
#> # ℹ 790 more rows
#> # ℹ 1 more variable: norm_coefficient <dbl>
```

<img src="man/figures/README-example_normalisation_plot-1.png"
style="width:100.0%" />

### Analysis

*analyse_TPP* fits sigmoidal melting curves are fitted to each unique
combination of protein, condition and replicate, and features of the
curve (including $T_m$) calculated. Significance statistics are
optionally calculated from the data:

- FDR-adjusted p-values are calculated for melting point differences
  from control conditions ($\Delta T_m$) for each
  replicate<sup>[1](#ref-Savitski_2014)</sup>

- Scaled F-scores and FDR-adjusted p-values are calculated from *NPARC*
  analysis for each protein<sup>[2](#ref-Childs_2019)</sup>

Default behaviour is to apply both analyses:

``` r
# Analyse twenty-protein data
twenty_prot_analysed <- 
  analyse_TPP(TPP_tbl = twenty_prot_normalised,
              control_name = "Control",
              p_value_methods = c("melting_point", "NPARC"))

# Analysed data and fitted curves
plot_melt(twenty_prot_analysed)
#> Warning: Removed 1000 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 3000 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
#> Warning: Removed 4000 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
```

<img src="man/figures/README-example_analysis-1.png" width="100%" />

All statistics are appended as new measurements to the *tibble*:

``` r
# Full table
twenty_prot_analysed
#> # A tibble: 800 × 25
#>    Protein_ID Condition Replicate  Temp Pep_N Match_N rel_quantity raw_quantity
#>    <chr>      <chr>     <chr>     <int> <dbl>   <dbl>        <dbl>        <dbl>
#>  1 Protein_A  Control   01           37    29      42       1          1401047.
#>  2 Protein_A  Control   01           47    29      42       1.15       1568061.
#>  3 Protein_A  Control   01           50    29      42       1.30       1889459 
#>  4 Protein_A  Control   01           59    29      42       0.0904      162054.
#>  5 Protein_A  Control   01           67    29      42       0.166       232101.
#>  6 Protein_A  Control   01           53    29      42       1.23       1801219.
#>  7 Protein_A  Control   01           56    29      42       0.450       630379.
#>  8 Protein_A  Control   01           44    29      42       1.03       1485400.
#>  9 Protein_A  Control   01           63    29      42       0.131       132437.
#> 10 Protein_A  Control   01           41    29      42       0.999      1464157.
#> # ℹ 790 more rows
#> # ℹ 17 more variables: norm_coefficient <dbl>, F_scaled <dbl>,
#> #   p_adj_NPARC <dbl>, a <dbl>, b <dbl>, plateau <dbl>, melt_point <dbl>,
#> #   infl_point <dbl>, slope <dbl>, R_sq <dbl>, Comparison <chr>,
#> #   diff_melt_point <dbl>, min_comparison_slope <dbl>, min_R_sq <dbl>,
#> #   min_slope <dbl>, max_control_plateau <dbl>, adj_pvalue <dbl>

# Statistics only
one_prot_stats <- 
  twenty_prot_analysed |>
  dplyr::filter(Protein_ID == "Protein_A", Condition != "Control") |>
  dplyr::select(Protein_ID, 
                Condition, 
                Replicate, 
                F_scaled, 
                p_adj_NPARC,
                melt_point,
                diff_melt_point,
                adj_pvalue
                ) |>
  dplyr::distinct()

# 1. Melting point
one_prot_stats |>
  dplyr::select(Protein_ID, 
                Condition, 
                Replicate,
                melt_point,
                diff_melt_point,
                adj_pvalue
                )
#> # A tibble: 2 × 6
#>   Protein_ID Condition Replicate melt_point diff_melt_point adj_pvalue
#>   <chr>      <chr>     <chr>          <dbl>           <dbl>      <dbl>
#> 1 Protein_A  Treated   01              56.0          0.0801      0.998
#> 2 Protein_A  Treated   02              56.0          0.0851      0.998

# 2. NPARC
one_prot_stats |>
  dplyr::select(Protein_ID, 
                Condition,
                F_scaled, 
                p_adj_NPARC,
                ) |>
  dplyr::distinct()
#> # A tibble: 1 × 4
#>   Protein_ID Condition F_scaled p_adj_NPARC
#>   <chr>      <chr>        <dbl>       <dbl>
#> 1 Protein_A  Treated      0.128       0.994
```

### Hit Identification

*get_TPP_hits* filters and summarises analysed data. Hits, by default,
have an NPARC-derived, FDR-adjusted p-value below
0.05,<sup>[2](#ref-Childs_2019)</sup> $\Delta T_m$ in the same direction
across replicates, and $\Delta T_m$ for comparisons against controls
greater than between controls.<sup>[1](#ref-Savitski_2014)</sup>

``` r
# Get hits from analysed twenty-protein data
twenty_prot_hits <- 
  get_TPP_hits(TPP_data = twenty_prot_analysed,
               hit_criteria = "default_hit_criteria",
               to_plot = TRUE,
               annotate = "melt_point")
#> --------------------
#> TPP Hit Identification
#> --------------------
#> Hit Criteria:
#>  NPARC_pvalue_threshold DTm_same_sign DTm_gt_Dcontrol
#>                    0.05          TRUE            TRUE
#> --------------------
#> Exporting hit data...
#> --------------------
#> TPP Export
#> --------------------
#> Export format: xlsx 
#> Saving TPP_hits.xlsx ...
#> Saved.
#> --------------------
#> Plotting hit melting curves...
#> 2 hits found.
#> --------------------

twenty_prot_hits
#> # A tibble: 2 × 10
#>   Protein_ID Condition Comparison         F_scaled p_adj_NPARC max_adj_pvalue
#>   <chr>      <chr>     <chr>                 <dbl>       <dbl>          <dbl>
#> 1 Protein_I  Treated   Treated_vs_Control     23.2   0.0000901      0.0000141
#> 2 Protein_C  Treated   Treated_vs_Control     14.1   0.00118        0.0981   
#> # ℹ 4 more variables: mean_melt_point <dbl>, mean_control_melt_point <dbl>,
#> #   mean_diff_melt_point <dbl>, mean_control_diff_melt_point <dbl>
```

<img src="man/figures/README-example_hits-1.png" style="width:100.0%" />

### Export

*get_TPP_hits* exported hit data by default, and full data can also be
exported as an excel .xlsx spreadsheet with *export_TPP*. Alternatives
are delimited text files (.csv, .tsv) or R data.

``` r
# Export as .xlsx
export_TPP(TPP_data = twenty_prot_analysed,
           file_name = "TPP_results.xlsx",
           format = "xlsx")
```

## References

<!-- ### Bibliography -->

<div id="refs" class="references csl-bib-body" entry-spacing="0"
line-spacing="2">

<div id="ref-Savitski_2014" class="csl-entry">

<span class="csl-left-margin">1.
</span><span class="csl-right-inline">Savitski, M. M. *et al.* [Tracking
cancer drugs in living cells by thermal profiling of the
proteome](https://doi.org/10.1126/science.1255784). *Science* **346**,
(2014).</span>

</div>

<div id="ref-Childs_2019" class="csl-entry">

<span class="csl-left-margin">2.
</span><span class="csl-right-inline">Childs, D. *et al.* [Nonparametric
analysis of thermal proteome profiles reveals novel drug-binding
proteins\*](https://doi.org/10.1074/mcp.tir119.001481). *Molecular &
Cellular Proteomics* **18**, 2506–2515 (2019).</span>

</div>

</div>
