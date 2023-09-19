
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
```

``` r
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
#> 1     t1       6.7056149       4.3714390    4.2409869      -4.3295817
#> 2    t10      -0.7187816      -0.7794280   -0.1848752      -0.7515719
#> 3   t100       0.9698969      -0.3729992   -0.4220562      -1.7705796
#> 4    t11      -8.1921507      -7.9036570   -5.4179821       3.6715388
#> 5    t12      -6.8151902      -6.2354807   -4.9236833       4.0297853
#> 6    t13       5.1227655       2.7839965   -0.7544528      -2.2194609
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1      -6.2570886   -2.7613480         0.65624630          0.2749146
#> 2      -0.8810548   -1.5957439         0.23747458          0.5849324
#> 3      -3.5535801   -2.1647978        -0.76145472         -0.4045134
#> 4       2.7472713   -0.7135728         0.79905161          0.4568037
#> 5       3.0099606    1.9836884        -0.04519136         -0.4181187
#> 6      -4.8975402   -2.5145352         0.39690857          0.3402024
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1      0.61856750        -1.30358161        -1.18267511    -0.730151517
#> 2      0.31561053        -0.50906818        -0.08611124    -0.473420949
#> 3     -0.09126971         0.08417716         0.55572570    -0.001290057
#> 4      0.55534859        -0.09578847         0.11519576     0.102301563
#> 5      0.26945938        -1.40919917        -0.75211615    -0.243613155
#> 6      0.86802357         0.33457711         0.58499789     0.676467689
```

We can use the function *plotData* to plot some simulated traits. Let’s
plot the phylo_G1 traits and environmental variable.

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

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

In this plot we can se how the simulated variances and covariances are
phylogenetically conserved. Let’s now plot the nonPhylo_G1 traits and
environmental variable.

    #> ℹ invalid tbl_tree object. Missing column: parent,node.
    #> ℹ invalid tbl_tree object. Missing column: parent,node.
    #> ℹ invalid tbl_tree object. Missing column: parent,node.
    #> ℹ invalid tbl_tree object. Missing column: parent,node.
    #> ℹ invalid tbl_tree object. Missing column: parent,node.
    #> ℹ invalid tbl_tree object. Missing column: parent,node.

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" /> We
see how these variables variances and covariances are not
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
    #> 1 nonPhylo_G1_trait1                 100             0.2806931
    #> 2 nonPhylo_G1_trait2                 100             0.1913426
    #> 3    phylo_G1_trait1                 100             0.9882768
    #> 4    phylo_G1_trait2                 100             0.9957512
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.719306946                             0
    #> 2               0.796306787                             0
    #> 3               0.009402907                             0
    #> 4              -0.002318771                             0
    #> [1] "Covariance results:"
    #>              trait_1            trait_2 number_observations total_correlation
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.89560368
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.14074178
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.06525730
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.08037902
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.04954784
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.83407987
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.16493454                  0.730669143
    #> 2              -0.14229545                  0.001553663
    #> 3              -0.05939504                 -0.005862255
    #> 4               0.09292234                 -0.012543316
    #> 5              -0.04802251                 -0.001525333
    #> 6               0.81906064                  0.015019230
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                0.00000000                       0.00000000
    #> 2                0.04444444                       0.04444444
    #> 3                0.48888889                       0.53333333
    #> 4                0.35555556                       0.35555556
    #> 5                0.62222222                       0.53333333
    #> 6                0.00000000                       0.00000000
    #>   p_value_non_phylogenetic_correlation
    #> 1                            0.0000000
    #> 2                            0.8888889
    #> 3                            0.4888889
    #> 4                            0.3555556
    #> 5                            0.8444444
    #> 6                            0.0000000

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2806931
#> 2 nonPhylo_G1_trait2                 100             0.1913426
#> 3    phylo_G1_trait1                 100             0.9882768
#> 4    phylo_G1_trait2                 100             0.9957512
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.719306946                             0
#> 2               0.796306787                             0
#> 3               0.009402907                             0
#> 4              -0.002318771                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.89560368
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.14074178
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.06525730
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.08037902
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.04954784
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.83407987
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.16493454                  0.730669143
#> 2              -0.14229545                  0.001553663
#> 3              -0.05939504                 -0.005862255
#> 4               0.09292234                 -0.012543316
#> 5              -0.04802251                 -0.001525333
#> 6               0.81906064                  0.015019230
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.00000000
#> 2                0.04444444                       0.04444444
#> 3                0.48888889                       0.53333333
#> 4                0.35555556                       0.35555556
#> 5                0.62222222                       0.53333333
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.8888889
#> 3                            0.4888889
#> 4                            0.3555556
#> 5                            0.8444444
#> 6                            0.0000000
```

We can now plot this results as trait networks using the
**plotNetworks**. Let’s first plot the phylogenetic part of the variance
and covariance.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  0                  0                  1                  1

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.82

Let’s now plot the non-phylogenetic variances and covariances.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  1                  1                  0                  0

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.73

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
    #> 1 nonPhylo_G1_trait1                 100             0.2229827
    #> 2 nonPhylo_G1_trait2                 100             0.2768939
    #> 3    phylo_G1_trait1                 100             0.9787892
    #> 4    phylo_G1_trait2                 100             0.9766271
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1              7.770173e-01                             0
    #> 2              7.113694e-01                             0
    #> 3              9.889198e-03                             0
    #> 4              3.786607e-05                             0
    #>   environmental_variable non_attributed_phylogenetic_variance
    #> 1           phylo_G1_env                            0.1415422
    #> 2           phylo_G1_env                            0.2327687
    #> 3           phylo_G1_env                            0.3735428
    #> 4           phylo_G1_env                            0.3178595
    #>   environmental_phylogenetic_variance labile_environmental_variance
    #> 1                          0.08144052                   0.044926500
    #> 2                          0.04412520                   0.034603480
    #> 3                          0.60524635                   0.006453788
    #> 4                          0.65876760                   0.007670498
    #>   residual_variance p_value_non_attributed_phylogenetic_variance
    #> 1        0.73209076                                            0
    #> 2        0.68850259                                            0
    #> 3        0.01475705                                            0
    #> 4        0.01570238                                            0
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
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.85941043
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.05096598
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.07153979
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.13592061
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.15413998
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.85699272
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.16462414                 0.6947862977
    #> 2               0.06127690                -0.0103109183
    #> 3              -0.07199339                 0.0004536005
    #> 4               0.14752057                -0.0115999560
    #> 5              -0.16176827                 0.0076282887
    #> 6               0.84294428                 0.0140484347
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                        0.0000000
    #> 2                 0.4888889                        0.4888889
    #> 3                 0.5777778                        0.6666667
    #> 4                 0.4444444                        0.4000000
    #> 5                 0.3555556                        0.3555556
    #> 6                 0.0000000                        0.0000000
    #>   p_value_non_phylogenetic_correlation environmental_variable
    #> 1                            0.0000000           phylo_G1_env
    #> 2                            0.5777778           phylo_G1_env
    #> 3                            0.9777778           phylo_G1_env
    #> 4                            0.5333333           phylo_G1_env
    #> 5                            0.8444444           phylo_G1_env
    #> 6                            0.0000000           phylo_G1_env
    #>   non_attributed_phylogenetic_correlation
    #> 1                              0.08208926
    #> 2                              0.06733262
    #> 3                             -0.18441722
    #> 4                              0.08448050
    #> 5                             -0.07971151
    #> 6                              0.21974229
    #>   environmental_phylogenetic_correlation labile_environmental_correlation
    #> 1                            0.082534876                     4.700253e-02
    #> 2                           -0.006055726                    -3.730804e-03
    #> 3                            0.112423824                    -1.013161e-02
    #> 4                            0.063040068                    -4.509113e-03
    #> 5                           -0.082056759                    -9.087947e-05
    #> 6                            0.623201996                     6.415768e-03
    #>   residual_correlation p_value_non_attributed_phylogenetic_correlation
    #> 1          0.647783770                                      0.00000000
    #> 2         -0.006580114                                      0.22222222
    #> 3          0.010585210                                      0.04444444
    #> 4         -0.007090842                                      0.57777778
    #> 5          0.007719168                                      0.40000000
    #> 6          0.007632667                                      0.00000000
    #>   p_value_environmental_phylogenetic_correlation
    #> 1                                      0.2222222
    #> 2                                      0.9777778
    #> 3                                      0.4888889
    #> 4                                      0.5333333
    #> 5                                      0.6222222
    #> 6                                      0.0000000
    #>   p_value_labile_environmental_correlation p_value_residual_correlation
    #> 1                               0.35555556                    0.0000000
    #> 2                               0.84444444                    0.5333333
    #> 3                               0.57777778                    0.5333333
    #> 4                               0.75555556                    0.7555556
    #> 5                               0.93333333                    0.8444444
    #> 6                               0.08888889                    0.0000000

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2229827
#> 2 nonPhylo_G1_trait2                 100             0.2768939
#> 3    phylo_G1_trait1                 100             0.9787892
#> 4    phylo_G1_trait2                 100             0.9766271
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1              7.770173e-01                             0
#> 2              7.113694e-01                             0
#> 3              9.889198e-03                             0
#> 4              3.786607e-05                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1415422
#> 2           phylo_G1_env                            0.2327687
#> 3           phylo_G1_env                            0.3735428
#> 4           phylo_G1_env                            0.3178595
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.08144052                   0.044926500
#> 2                          0.04412520                   0.034603480
#> 3                          0.60524635                   0.006453788
#> 4                          0.65876760                   0.007670498
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.73209076                                            0
#> 2        0.68850259                                            0
#> 3        0.01475705                                            0
#> 4        0.01570238                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.85941043
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.05096598
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.07153979
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        0.13592061
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.15413998
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.85699272
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.16462414                 0.6947862977
#> 2               0.06127690                -0.0103109183
#> 3              -0.07199339                 0.0004536005
#> 4               0.14752057                -0.0115999560
#> 5              -0.16176827                 0.0076282887
#> 6               0.84294428                 0.0140484347
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                        0.0000000
#> 2                 0.4888889                        0.4888889
#> 3                 0.5777778                        0.6666667
#> 4                 0.4444444                        0.4000000
#> 5                 0.3555556                        0.3555556
#> 6                 0.0000000                        0.0000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                            0.0000000           phylo_G1_env
#> 2                            0.5777778           phylo_G1_env
#> 3                            0.9777778           phylo_G1_env
#> 4                            0.5333333           phylo_G1_env
#> 5                            0.8444444           phylo_G1_env
#> 6                            0.0000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.08208926
#> 2                              0.06733262
#> 3                             -0.18441722
#> 4                              0.08448050
#> 5                             -0.07971151
#> 6                              0.21974229
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.082534876                     4.700253e-02
#> 2                           -0.006055726                    -3.730804e-03
#> 3                            0.112423824                    -1.013161e-02
#> 4                            0.063040068                    -4.509113e-03
#> 5                           -0.082056759                    -9.087947e-05
#> 6                            0.623201996                     6.415768e-03
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.647783770                                      0.00000000
#> 2         -0.006580114                                      0.22222222
#> 3          0.010585210                                      0.04444444
#> 4         -0.007090842                                      0.57777778
#> 5          0.007719168                                      0.40000000
#> 6          0.007632667                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.2222222
#> 2                                      0.9777778
#> 3                                      0.4888889
#> 4                                      0.5333333
#> 5                                      0.6222222
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                               0.35555556                    0.0000000
#> 2                               0.84444444                    0.5333333
#> 3                               0.57777778                    0.5333333
#> 4                               0.75555556                    0.7555556
#> 5                               0.93333333                    0.8444444
#> 6                               0.08888889                    0.0000000
```

Let’s plot the phylogenetic environmental variances and covariances as
an example.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  0                  0                  1                  1

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.62

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
    #> 1    phylo_G1_env                  80             0.9921964
    #> 2 phylo_G1_trait1                  80             0.9764975
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.007803577                             0
    #> 2               0.012821580                             0
    #> [1] "Covariance results:"
    #>        trait_1         trait_2 number_observations total_correlation
    #> 1 phylo_G1_env phylo_G1_trait1                  80         0.7541053
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1                0.7432015                   0.01090379
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
    #> 1 phylo_G1_trait1 100    80   20 1.123189
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.130978
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.147844
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.157026
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.159605
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.110604
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.164395
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.109199
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.100836
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.140279
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.167397
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.116503
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.116675
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.182371
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.120062
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.145348
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.092878
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.157411
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.135907
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.111319
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.111213
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.153975
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.10573
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.11178
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.185329
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA  NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.1466
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.120348
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.114945
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.09509
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.124258
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1       6.7056149
#> 2    t10      -0.3547771
#> 3   t100       0.9698969
#> 4    t11      -8.1921507
#> 5    t12      -6.8151902
#> 6    t13       5.1227655
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
#>  -8.192  -1.780   1.734   1.033   3.382   7.124      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -8.1922 -1.7804  1.6746  0.9553  3.2451  7.1244
```
