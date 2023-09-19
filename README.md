
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TrEvol

<!-- badges: start -->
<!-- badges: end -->

The goal of TrEvol is to provide easy to use functions to study
functional trait syndromes evolutionary patterns. More concretly, this
package is focused on separating trait variance and the covariance
(represented as correlation) among pairs of traits into a phylogenetic
and a non-phylogenetic component. This allows to quantify the amount of
phylogenetic conservatism (phylogenetic signal) and evolutionary
lability in traits and their correlation.

Moreover, the package allows to consider one environmental variable when
doing that, which allows to further separate trait variances and
correlation among pairs of traits into four components, quantifying its
relationship with the phylogeny and an environmental variable of
interest. The first component is the *non-attributed phylogenetic
variances and correlation*, and represents the variance and correlation
that is related only to the phylogeny (not related to the environmental
variable). The second component is the *environmental phylogenetic
variances and correlation*, representing the variances and correlation
related both to the phylogeny and an environmental variable. The third
component is the *labile environmental variances and correlation*, and
represents the amount of variance and correlation that is only related
to the environmental variable (not related to the phylogeny). Finally,
the fourth component is the *residual variances and correlation*, which
is not related to the phylogeny nor to an environmental variable of
interest. Total correlation among traits is also calculated.

The package allows then to characterize variance and covariance patterns
for a group of traits. To represent these results, the package use trait
networks, showing the variance components in the nodes and the
correlation among traits as edges. The package allows to display network
metrics describing the structure of the network. Networks can be
represented for each of the variance and covariance components
previously described, and they can be compared.

The package also includes an imputation framework which uses random
forest to predict missing values including phylogenetic and
environmental information as well as the relationship between traits.
This imputation framework optimizes the use of the information by
selecting the elements that are expected to be bette predictors and uses
them to impute missing values in a given dataset.

## Installation

You can install the development version of TrEvol from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pablosanchezmart/TrEvol")
```

## Example

### Simulate example data

This is a basic example showing how to use the package using simulated
data. First of all, let’s simulate some data to work with. To do so, we
will use the **simulateDataSet** function from the package using the
default parameters.

``` r
library(TrEvol)

simulated.data <- simulateDataSet(number_observations = 100)
```

The **simulateDataSet** function simulates a set of correlated variables
and a phylogeny. The variance-covariance matrix used to simulate traits
variances (diagonal) and correlations (off-diagonal) is reported by the
function and can be modified using the *vcv_matrix* argument.

``` r
simulated.data$vcv_matrix
#>           G1_trait1 G1_trait2 G1_env G2_trait1 G2_trait2 G2_env
#> G1_trait1       1.0       0.9    0.8       0.0       0.1    0.2
#> G1_trait2       0.9       1.0    0.8       0.0       0.1    0.2
#> G1_env          0.8       0.8    1.0       0.0       0.1    0.2
#> G2_trait1       0.0       0.0    0.0       1.0       0.9    0.8
#> G2_trait2       0.1       0.1    0.1       0.9       1.0    0.8
#> G2_env          0.2       0.2    0.2       0.8       0.8    1.0
```

The **simulateDataSet** function produces two sets of variables. The
firs set uses the variance-covariance matrix to simulate
phylogenetically structured data under a Brownian motion model of
evolution. These variables are preceded by “phylo\_”. To do so, the
function simulates a phylogeny with the same number of tips as number of
observations to be simulated (set in the number_observations argument).
A phylogeny can be introduced manually using the *phylogeny* argument.
The second set of variables are simulated in a non-phylogenetically
structured manner. So their variances and covariances will not present a
phylogenetic component. These variables are preceded by “nonPhylo\_”.
Each of the sets have 6 variables, 4 will be considered traits and one
will be considered as an enviornmental variable. Names of each variable
show the expected correlation group (G1 or G2, for each set of traits).
So we expect phylo G1 to be correlated among them and to present
variances and covariances related to the phylogeny and to the
phylo_G1_env environmental variable.

``` r
head(simulated.data$data)
#>   animal phylo_G1_trait1 phylo_G1_trait2 phylo_G1_env phylo_G2_trait1
#> 1     t1      -0.7791497       -1.146632   -5.3511190       2.3733229
#> 2    t10       4.0229772        5.196852    4.4245836      -4.8983153
#> 3   t100       2.2433712        2.452164    3.1806707       3.2117123
#> 4    t11       2.6256106        4.321767    2.6815642      -0.3606371
#> 5    t12       0.2905713       -1.513021   -3.4163089      -0.1829847
#> 6    t13       1.3612768        2.599177    0.8945128      -0.8599589
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1       1.7160880   0.51430748        -0.06109192        -0.07606438
#> 2      -1.6928676  -3.80130077        -0.33756837        -0.37502361
#> 3       3.2529266   1.09927359         0.76522664         0.71141334
#> 4       2.3203191   0.74755243         0.05722483         0.98416368
#> 5       0.5091803  -2.57958370        -1.23241135        -1.36138756
#> 6       0.4986816  -0.05204808         0.47823186        -0.19599444
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1     -0.74035540         1.78209639          1.3993626       1.1472065
#> 2      0.40886206        -1.55866634         -1.0399105      -0.5968838
#> 3      0.20589097         0.36488928          0.4729491       0.8516450
#> 4      0.02203621        -0.76450168         -1.4220850      -1.7166284
#> 5     -0.67643787         0.74078276          0.5775648       1.0457082
#> 6     -0.29814243        -0.05313711          0.1790064      -0.6197324
```

We can use the function *plotData* to plot some simulated traits. Let’s
plot the phylo_G1 traits and environmental variable.

``` r
plotData(variables = c("phylo_G1_trait1", "phylo_G1_trait2", "phylo_G1_env"),
         dataset = simulated.data$data,
         phylogeny = simulated.data$phylogeny,
         terminal.taxon = "animal" 
          )
#> Registered S3 methods overwritten by 'treeio':
#>   method              from    
#>   MRCA.phylo          tidytree
#>   MRCA.treedata       tidytree
#>   Nnode.treedata      tidytree
#>   Ntip.treedata       tidytree
#>   ancestor.phylo      tidytree
#>   ancestor.treedata   tidytree
#>   child.phylo         tidytree
#>   child.treedata      tidytree
#>   full_join.phylo     tidytree
#>   full_join.treedata  tidytree
#>   groupClade.phylo    tidytree
#>   groupClade.treedata tidytree
#>   groupOTU.phylo      tidytree
#>   groupOTU.treedata   tidytree
#>   is.rooted.treedata  tidytree
#>   nodeid.phylo        tidytree
#>   nodeid.treedata     tidytree
#>   nodelab.phylo       tidytree
#>   nodelab.treedata    tidytree
#>   offspring.phylo     tidytree
#>   offspring.treedata  tidytree
#>   parent.phylo        tidytree
#>   parent.treedata     tidytree
#>   root.treedata       tidytree
#>   rootnode.phylo      tidytree
#>   sibling.phylo       tidytree
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-4-2.png" width="100%" />

In this plot we can se how the simulated variances and covariances are
phylogenetically conserved. Let’s now plot the nonPhylo_G1 traits and
environmental variable.

``` r
plotData(variables = c("nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "nonPhylo_G1_env"),
         dataset = simulated.data$data,
         phylogeny = simulated.data$phylogeny,
         terminal.taxon = "animal" 
          )
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
#> ℹ invalid tbl_tree object. Missing column: parent,node.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />
We see how these variables variances and covariances are not
phylogenetically conserved.

### Compute variances and covariances for simulated data

Now we can use the function **computeVarianceCovariancePartition** to
report the correlation structure between traits.

Let’s imagine that we are only interested in computing the phylogenetic
and non-phylogenetic variances and covariances. So, for now, let’s nos
include an environmental variable in the function. In this case,
phylogenetic variances and correlation and non-phylogenetic variances
and correlation will be reported, jointly with total correlation. We
will use the default model specifications, but users can modify this by
using the defineModelsSpecifications function to create an object that
then can be introduced in the *model_specifications* argument of the
**computeVarianceCovariancePartition** function. This will be specially
needed with complex data, where number of iterations, burning and
thinning will need to be increased. Users can also explore using
different priors.

``` r
variance_covariance_results.list <- computeVarianceCovariancePartition(
  traits =  c("nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "phylo_G1_trait1", "phylo_G1_trait2"),
  dataset = simulated.data$data,
  terminal_taxa = "animal",
  phylogeny = simulated.data$phylogeny
  )
#> [1] "Running variance-covariance calculation: bi_nonPhylo_G1_trait1_nonPhylo_G1_trait2"
#> [1] "Using default model specificatios. Use defineModelsSpecifications() output on model_specifications argument to set them manually."
#> [1] "Running variance-covariance calculation: bi_nonPhylo_G1_trait1_phylo_G1_trait1"
#> [1] "Running variance-covariance calculation: bi_nonPhylo_G1_trait1_phylo_G1_trait2"
#> [1] "Running variance-covariance calculation: bi_nonPhylo_G1_trait2_phylo_G1_trait1"
#> [1] "Running variance-covariance calculation: bi_nonPhylo_G1_trait2_phylo_G1_trait2"
#> [1] "Running variance-covariance calculation: bi_phylo_G1_trait1_phylo_G1_trait2"
#> [1] "Model structure used:"
#>                                       type
#> 1 bi_nonPhylo_G1_trait1_nonPhylo_G1_trait2
#> 2    bi_nonPhylo_G1_trait1_phylo_G1_trait1
#> 3    bi_nonPhylo_G1_trait1_phylo_G1_trait2
#> 4    bi_nonPhylo_G1_trait2_phylo_G1_trait1
#> 5    bi_nonPhylo_G1_trait2_phylo_G1_trait2
#> 6       bi_phylo_G1_trait1_phylo_G1_trait2
#>                                   traits             trait1             trait2
#> 1 nonPhylo_G1_trait1, nonPhylo_G1_trait2 nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 2    nonPhylo_G1_trait1, phylo_G1_trait1 nonPhylo_G1_trait1    phylo_G1_trait1
#> 3    nonPhylo_G1_trait1, phylo_G1_trait2 nonPhylo_G1_trait1    phylo_G1_trait2
#> 4    nonPhylo_G1_trait2, phylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1
#> 5    nonPhylo_G1_trait2, phylo_G1_trait2 nonPhylo_G1_trait2    phylo_G1_trait2
#> 6       phylo_G1_trait1, phylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2
#>                                                  fix.frml           ran.frml
#> 1 cbind(nonPhylo_G1_trait1, nonPhylo_G1_trait2) ~ trait-1 ~ us(trait):animal
#> 2    cbind(nonPhylo_G1_trait1, phylo_G1_trait1) ~ trait-1 ~ us(trait):animal
#> 3    cbind(nonPhylo_G1_trait1, phylo_G1_trait2) ~ trait-1 ~ us(trait):animal
#> 4    cbind(nonPhylo_G1_trait2, phylo_G1_trait1) ~ trait-1 ~ us(trait):animal
#> 5    cbind(nonPhylo_G1_trait2, phylo_G1_trait2) ~ trait-1 ~ us(trait):animal
#> 6       cbind(phylo_G1_trait1, phylo_G1_trait2) ~ trait-1 ~ us(trait):animal
#> [1] "Variance results:"
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.1459568
#> 2 nonPhylo_G1_trait2                 100             0.1945975
#> 3    phylo_G1_trait1                 100             0.9762921
#> 4    phylo_G1_trait2                 100             0.9705034
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.854043218                             0
#> 2               0.806898919                             0
#> 3              -0.000827434                             0
#> 4              -0.027640206                             0
#> [1] "Covariance results:"
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.83598311
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.06711897
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.12929594
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.15676110
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.06674714
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90044476
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1              0.115135236                   0.72084787
#> 2              0.006262329                  -0.07338130
#> 3              0.197268684                  -0.06797275
#> 4              0.223944375                  -0.06718327
#> 5              0.004550602                  -0.07129774
#> 6              0.880192803                   0.02025196
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                       0.00000000
#> 2                 0.3111111                       0.75555556
#> 3                 0.3555556                       0.08888889
#> 4                 0.0000000                       0.00000000
#> 5                 0.2666667                       0.75555556
#> 6                 0.0000000                       0.00000000
#>   p_value_non_phylogenetic_correlation
#> 1                           0.00000000
#> 2                           0.00000000
#> 3                           0.04444444
#> 4                           0.00000000
#> 5                           0.00000000
#> 6                           0.00000000
```

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.1459568
#> 2 nonPhylo_G1_trait2                 100             0.1945975
#> 3    phylo_G1_trait1                 100             0.9762921
#> 4    phylo_G1_trait2                 100             0.9705034
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.854043218                             0
#> 2               0.806898919                             0
#> 3              -0.000827434                             0
#> 4              -0.027640206                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.83598311
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.06711897
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.12929594
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.15676110
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.06674714
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90044476
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1              0.115135236                   0.72084787
#> 2              0.006262329                  -0.07338130
#> 3              0.197268684                  -0.06797275
#> 4              0.223944375                  -0.06718327
#> 5              0.004550602                  -0.07129774
#> 6              0.880192803                   0.02025196
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                       0.00000000
#> 2                 0.3111111                       0.75555556
#> 3                 0.3555556                       0.08888889
#> 4                 0.0000000                       0.00000000
#> 5                 0.2666667                       0.75555556
#> 6                 0.0000000                       0.00000000
#>   p_value_non_phylogenetic_correlation
#> 1                           0.00000000
#> 2                           0.00000000
#> 3                           0.04444444
#> 4                           0.00000000
#> 5                           0.00000000
#> 6                           0.00000000
```

We can now plot this results as trait networks using the
**plotNetworks**. Let’s first plot the phylogenetic part of the variance
and covariance.

``` r
plotNetwork(variance_results = variance_covariance_results.list$varianceResults,
            correlation_results = variance_covariance_results.list$covarianceResults,
            variance_type = "phylogenetic_variance", 
            correlation_type = "phylogenetic_correlation", 
            edge_label = T, 
            only_significant = T, 
            correlation_threshold = 0.3)
#> [1] "Node degree:"
#> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
#>                  0                  0                  1                  1
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.88

Let’s now plot the non-phylogenetic variances and covariances.

``` r
plotNetwork(variance_results = variance_covariance_results.list$varianceResults,
            correlation_results = variance_covariance_results.list$covarianceResults,
            variance_type = "non_phylogenetic_variance", 
            correlation_type = "non_phylogenetic_correlation", 
            edge_label = T, 
            only_significant = T, 
            correlation_threshold = 0.3)
#> [1] "Node degree:"
#> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
#>                  1                  1                  0                  0
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.72

We observe how much of the variance and covariance in the phylo_G1
traits is related to the phylogeny, and that the variances and
covariances in nonPhylo_G1 traits are not related to the phylogeny, as
expected.

### Compute variances and covariances for simulated data including an enviornmental variable

Let’s now ilustrate how to include an environmental variable. This will
allow us to calculate four different variance and covariance components:
the *non-attributed phylogenetic variance and covariance*, the
*phylogenetic environmental variance and covariance*, the *labile
environmental variance and covariance* and the *residual variance and
covariance*. Let’s use the phylo_G1_env variable as an enviornmental
variable, which is expected to be related to variances and covariances
of phylo_G1 traits in a phylogenetically conserved manner (i.e.,
phylogenetic environmental variance and covariance).

``` r
variance_covariance_environment_results.list <- computeVarianceCovariancePartition(
  traits =  c("nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "phylo_G1_trait1", "phylo_G1_trait2"),
  environmental_variable = "phylo_G1_env",
  dataset = simulated.data$data,
  terminal_taxa = "animal",
  phylogeny = simulated.data$phylogeny
  )
#> [1] "Running variance-covariance calculation: tri_nonPhylo_G1_trait1_nonPhylo_G1_trait2_phylo_G1_env"
#> [1] "Using default model specificatios. Use defineModelsSpecifications() output on model_specifications argument to set them manually."
#> [1] "Running variance-covariance calculation: tri_nonPhylo_G1_trait1_phylo_G1_trait1_phylo_G1_env"
#> [1] "Running variance-covariance calculation: tri_nonPhylo_G1_trait1_phylo_G1_trait2_phylo_G1_env"
#> [1] "Running variance-covariance calculation: tri_nonPhylo_G1_trait2_phylo_G1_trait1_phylo_G1_env"
#> [1] "Running variance-covariance calculation: tri_nonPhylo_G1_trait2_phylo_G1_trait2_phylo_G1_env"
#> [1] "Running variance-covariance calculation: tri_phylo_G1_trait1_phylo_G1_trait2_phylo_G1_env"
#> [1] "Model structure used:"
#>                                                     type
#> 1 tri_nonPhylo_G1_trait1_nonPhylo_G1_trait2_phylo_G1_env
#> 2    tri_nonPhylo_G1_trait1_phylo_G1_trait1_phylo_G1_env
#> 3    tri_nonPhylo_G1_trait1_phylo_G1_trait2_phylo_G1_env
#> 4    tri_nonPhylo_G1_trait2_phylo_G1_trait1_phylo_G1_env
#> 5    tri_nonPhylo_G1_trait2_phylo_G1_trait2_phylo_G1_env
#> 6       tri_phylo_G1_trait1_phylo_G1_trait2_phylo_G1_env
#>                                   traits             trait1             trait2
#> 1 nonPhylo_G1_trait1, nonPhylo_G1_trait2 nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 2    nonPhylo_G1_trait1, phylo_G1_trait1 nonPhylo_G1_trait1    phylo_G1_trait1
#> 3    nonPhylo_G1_trait1, phylo_G1_trait2 nonPhylo_G1_trait1    phylo_G1_trait2
#> 4    nonPhylo_G1_trait2, phylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1
#> 5    nonPhylo_G1_trait2, phylo_G1_trait2 nonPhylo_G1_trait2    phylo_G1_trait2
#> 6       phylo_G1_trait1, phylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2
#>                                                                fix.frml
#> 1 cbind(nonPhylo_G1_trait1, nonPhylo_G1_trait2, phylo_G1_env) ~ trait-1
#> 2    cbind(nonPhylo_G1_trait1, phylo_G1_trait1, phylo_G1_env) ~ trait-1
#> 3    cbind(nonPhylo_G1_trait1, phylo_G1_trait2, phylo_G1_env) ~ trait-1
#> 4    cbind(nonPhylo_G1_trait2, phylo_G1_trait1, phylo_G1_env) ~ trait-1
#> 5    cbind(nonPhylo_G1_trait2, phylo_G1_trait2, phylo_G1_env) ~ trait-1
#> 6       cbind(phylo_G1_trait1, phylo_G1_trait2, phylo_G1_env) ~ trait-1
#>             ran.frml
#> 1 ~ us(trait):animal
#> 2 ~ us(trait):animal
#> 3 ~ us(trait):animal
#> 4 ~ us(trait):animal
#> 5 ~ us(trait):animal
#> 6 ~ us(trait):animal
#> [1] "Variance results:"
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2148143
#> 2 nonPhylo_G1_trait2                 100             0.3289702
#> 3    phylo_G1_trait1                 100             0.9717023
#> 4    phylo_G1_trait2                 100             0.9691645
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.785185743                             0
#> 2               0.673843568                             0
#> 3               0.004717517                             0
#> 4              -0.006709929                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1904739
#> 2           phylo_G1_env                            0.3101555
#> 3           phylo_G1_env                            0.3680740
#> 4           phylo_G1_env                            0.3622588
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.02434037                    0.19230636
#> 2                          0.01881471                    0.11879175
#> 3                          0.60362826                    0.01177801
#> 4                          0.60690578                    0.01672543
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.59287938                                            0
#> 2        0.55223805                                            0
#> 3        0.01651973                                            0
#> 4        0.01411004                                            0
#>   p_value_environmental_phylogenetic_variance
#> 1                                           0
#> 2                                           0
#> 3                                           0
#> 4                                           0
#>   p_value_labile_environmental_variance
#> 1                                     0
#> 2                                     0
#> 3                                     0
#> 4                                     0
#> [1] "Covariance results:"
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.80917303
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.14230027
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.01283962
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.15123034
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        0.03540373
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.92086136
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.12834436                   0.68082867
#> 2               0.19154294                  -0.04924266
#> 3               0.03240365                  -0.01956403
#> 4               0.21803274                  -0.06680240
#> 5               0.09977065                  -0.06436693
#> 6               0.90091781                   0.01994355
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.04444444
#> 2                0.48888889                       0.17777778
#> 3                0.93333333                       0.80000000
#> 4                0.04444444                       0.00000000
#> 5                0.71111111                       0.35555556
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                           0.00000000           phylo_G1_env
#> 2                           0.04444444           phylo_G1_env
#> 3                           0.44444444           phylo_G1_env
#> 4                           0.00000000           phylo_G1_env
#> 5                           0.08888889           phylo_G1_env
#> 6                           0.00000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                               0.1107952
#> 2                               0.1618155
#> 3                               0.0348201
#> 4                               0.2122651
#> 5                               0.1379126
#> 6                               0.3030154
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.017549173                      0.168845352
#> 2                            0.029727421                     -0.023568489
#> 3                           -0.002416448                     -0.009596926
#> 4                            0.005767639                     -0.045304240
#> 5                           -0.038141923                     -0.045258116
#> 6                            0.597902401                      0.013978816
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.511983316                                      0.04444444
#> 2         -0.025674175                                      0.04444444
#> 3         -0.009967103                                      0.53333333
#> 4         -0.021498163                                      0.00000000
#> 5         -0.019108810                                      0.00000000
#> 6          0.005964734                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.5777778
#> 2                                      0.9333333
#> 3                                      0.9777778
#> 4                                      0.9777778
#> 5                                      0.5333333
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                               0.00000000                   0.00000000
#> 2                               0.13333333                   0.08888889
#> 3                               0.53333333                   0.57777778
#> 4                               0.22222222                   0.44444444
#> 5                               0.08888889                   0.31111111
#> 6                               0.04444444                   0.04444444
```

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2148143
#> 2 nonPhylo_G1_trait2                 100             0.3289702
#> 3    phylo_G1_trait1                 100             0.9717023
#> 4    phylo_G1_trait2                 100             0.9691645
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.785185743                             0
#> 2               0.673843568                             0
#> 3               0.004717517                             0
#> 4              -0.006709929                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1904739
#> 2           phylo_G1_env                            0.3101555
#> 3           phylo_G1_env                            0.3680740
#> 4           phylo_G1_env                            0.3622588
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.02434037                    0.19230636
#> 2                          0.01881471                    0.11879175
#> 3                          0.60362826                    0.01177801
#> 4                          0.60690578                    0.01672543
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.59287938                                            0
#> 2        0.55223805                                            0
#> 3        0.01651973                                            0
#> 4        0.01411004                                            0
#>   p_value_environmental_phylogenetic_variance
#> 1                                           0
#> 2                                           0
#> 3                                           0
#> 4                                           0
#>   p_value_labile_environmental_variance
#> 1                                     0
#> 2                                     0
#> 3                                     0
#> 4                                     0
```

``` r
variance_covariance_environment_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.80917303
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.14230027
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.01283962
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.15123034
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        0.03540373
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.92086136
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.12834436                   0.68082867
#> 2               0.19154294                  -0.04924266
#> 3               0.03240365                  -0.01956403
#> 4               0.21803274                  -0.06680240
#> 5               0.09977065                  -0.06436693
#> 6               0.90091781                   0.01994355
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.04444444
#> 2                0.48888889                       0.17777778
#> 3                0.93333333                       0.80000000
#> 4                0.04444444                       0.00000000
#> 5                0.71111111                       0.35555556
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                           0.00000000           phylo_G1_env
#> 2                           0.04444444           phylo_G1_env
#> 3                           0.44444444           phylo_G1_env
#> 4                           0.00000000           phylo_G1_env
#> 5                           0.08888889           phylo_G1_env
#> 6                           0.00000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                               0.1107952
#> 2                               0.1618155
#> 3                               0.0348201
#> 4                               0.2122651
#> 5                               0.1379126
#> 6                               0.3030154
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.017549173                      0.168845352
#> 2                            0.029727421                     -0.023568489
#> 3                           -0.002416448                     -0.009596926
#> 4                            0.005767639                     -0.045304240
#> 5                           -0.038141923                     -0.045258116
#> 6                            0.597902401                      0.013978816
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.511983316                                      0.04444444
#> 2         -0.025674175                                      0.04444444
#> 3         -0.009967103                                      0.53333333
#> 4         -0.021498163                                      0.00000000
#> 5         -0.019108810                                      0.00000000
#> 6          0.005964734                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.5777778
#> 2                                      0.9333333
#> 3                                      0.9777778
#> 4                                      0.9777778
#> 5                                      0.5333333
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                               0.00000000                   0.00000000
#> 2                               0.13333333                   0.08888889
#> 3                               0.53333333                   0.57777778
#> 4                               0.22222222                   0.44444444
#> 5                               0.08888889                   0.31111111
#> 6                               0.04444444                   0.04444444
```

Let’s plot the phylogenetic environmental variances and covariances as
an example.

``` r
plotNetwork(variance_results = variance_covariance_environment_results.list$varianceResults,
            correlation_results = variance_covariance_environment_results.list$covarianceResults,
            variance_type = "environmental_phylogenetic_variance", 
            correlation_type = "environmental_phylogenetic_correlation", 
            edge_label = T, 
            only_significant = T, 
            correlation_threshold = 0.3)
#> [1] "Node degree:"
#> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
#>                  0                  0                  1                  1
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.6

We can see how part of the pattern of variances and covariances of the
traits phylo_G1 are attributed to the phylogenetically conserved effect
of the environmental variable phylo_G1_env.

### Impute missing values

Now, let’s use this information to inform an imputation framework
performed by the function **imputeTraits** of this package. As our
simulated data does not have NAs, let’s first produce some NAs (20%) in
a variable to impute. Then, we will run the imputeTraits to perform
predictions using the variance-covariance structures reported before to
inform the process.

``` r
to_impute.data <- simulated.data$data
to_impute.data[, "phylo_G1_trait1"] <- missForest::prodNA(as.data.frame(simulated.data$data[, "phylo_G1_trait1"]), 0.2)

imputed.data <- imputeTraits(variables_to_impute = "phylo_G1_trait1",
                             dataset = to_impute.data, 
                             terminal_taxon = "animal", 
                             phylogeny = simulated.data$phylogeny,
                             predictors = "phylo_G1_env")
#> [1] "Running variance-covariance calculation: bi_phylo_G1_env_phylo_G1_trait1"
#> [1] "Model structure used:"
#>                              type                        traits       trait1
#> 1 bi_phylo_G1_env_phylo_G1_trait1 phylo_G1_env, phylo_G1_trait1 phylo_G1_env
#>            trait2                                       fix.frml
#> 1 phylo_G1_trait1 cbind(phylo_G1_env, phylo_G1_trait1) ~ trait-1
#>             ran.frml
#> 1 ~ us(trait):animal
#> [1] "Variance results:"
#>             trait number_observations phylogenetic_variance
#> 1    phylo_G1_env                  80             0.9565321
#> 2 phylo_G1_trait1                  80             0.9551053
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1                0.04346794                             0
#> 2                0.04793374                             0
#> [1] "Covariance results:"
#>        trait_1         trait_2 number_observations total_correlation
#> 1 phylo_G1_env phylo_G1_trait1                  80         0.8116482
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1                0.7694984                   0.04214977
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                         0                                0
#>   p_value_non_phylogenetic_correlation
#> 1                                    0
#> Loading required package: missForest
#> Loading required package: doRNG
#> Loading required package: foreach
#> Loading required package: rngtools
#> Loading required package: randomForest
#> randomForest 4.7-1.1
#> Type rfNews() to see new features/changes/bug fixes.
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.626723
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.647177
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.645399
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.646647
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.647881
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.607152
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.624049
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.666394
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.678829
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.655385
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.644359
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.64801
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.65064
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.68197
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.59831
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.646471
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.629096
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.659916
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.645022
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.669383
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.646394
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.639522
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.642611
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.68058
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.676642
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.612633
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.621294
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.65838
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.674563
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.618299
#>                                                                                                                                                                                                                            Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, phylo_G1_env

head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1      -0.7791497
#> 2    t10       4.5424677
#> 3   t100       2.4088649
#> 4    t11       2.6256106
#> 5    t12      -1.4174879
#> 6    t13       1.3612768
```

Let’s assign the imputed data as a new column in the dataset.

``` r
imputed.data <- merge(to_impute.data, imputed.data$round3$ximp, 
                      suffixes = c("", "_imputed"), 
                      by = "animal", 
                      all.x = T)
```

Let’s look at the summary of the imputed data and compare it with the
complete data.

Before imputation:

``` r
summary(imputed.data$phylo_G1_trait1)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#>  -4.545   1.005   2.876   2.488   4.175   6.340      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1imputed)
#> Length  Class   Mode 
#>      0   NULL   NULL
```
