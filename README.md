
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

## Simulate example data

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
#> 1     t1        1.843047       1.1122000    4.0087633      -0.2670343
#> 2    t10       -1.052142      -0.9919674   -3.3640685      -1.9877084
#> 3   t100       -1.682892      -1.5868368    0.8166605       1.5585355
#> 4    t11        3.447466       4.5056798    1.6619987      -3.5886639
#> 5    t12       -2.662426      -3.0764030   -1.6034489       2.9807050
#> 6    t13       -6.524720      -4.3946440   -1.3897919       0.3552793
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1      -0.2139250  -1.25462311         0.04667906         -0.1486699
#> 2      -1.3211920  -0.04512673         0.31725221          0.1768965
#> 3       0.4378832   2.31746877         0.06222407          0.4761982
#> 4      -2.0696231  -2.71727580        -1.83314082         -1.7253070
#> 5       1.0437484   1.17455735        -0.80000481         -1.5477175
#> 6      -1.0721860   0.43346537         0.52935238          0.7040617
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1      -0.3571845         -0.8963399         -1.3762793     -0.09617784
#> 2      -0.5740196         -1.2752680         -1.0810088     -1.91180864
#> 3      -0.6457899         -0.9417485         -0.7276413     -0.82339751
#> 4      -1.1420439          1.0109764          0.4963721     -0.10093743
#> 5      -0.5959015         -2.3025604         -2.5555673     -0.03558588
#> 6       0.1622648          2.1097471          2.1463683      2.15046869
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

## Compute variances and covariances for simulated data

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
    #> 1 nonPhylo_G1_trait1                 100             0.3645125
    #> 2 nonPhylo_G1_trait2                 100             0.3671487
    #> 3    phylo_G1_trait1                 100             0.9899829
    #> 4    phylo_G1_trait2                 100             0.9869998
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.635487516                             0
    #> 2               0.626416258                             0
    #> 3               0.007813634                             0
    #> 4               0.004666303                             0
    #> [1] "Covariance results:"
    #>              trait_1            trait_2 number_observations total_correlation
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.88790784
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.03648434
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.10030051
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.30032143
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.31187340
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90849496
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.38660182                 0.5013060193
    #> 2               0.02618614                 0.0102981933
    #> 3              -0.10998772                 0.0096872090
    #> 4              -0.29960088                -0.0007205544
    #> 5              -0.33078317                 0.0189097688
    #> 6               0.89746122                 0.0110337423
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                       0.00000000
    #> 2                 0.6666667                       0.84444444
    #> 3                 0.5333333                       0.48888889
    #> 4                 0.0000000                       0.04444444
    #> 5                 0.0000000                       0.00000000
    #> 6                 0.0000000                       0.00000000
    #>   p_value_non_phylogenetic_correlation
    #> 1                            0.0000000
    #> 2                            0.4000000
    #> 3                            0.5777778
    #> 4                            0.6222222
    #> 5                            0.4000000
    #> 6                            0.0000000

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.3645125
#> 2 nonPhylo_G1_trait2                 100             0.3671487
#> 3    phylo_G1_trait1                 100             0.9899829
#> 4    phylo_G1_trait2                 100             0.9869998
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.635487516                             0
#> 2               0.626416258                             0
#> 3               0.007813634                             0
#> 4               0.004666303                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.88790784
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.03648434
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.10030051
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.30032143
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.31187340
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90849496
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.38660182                 0.5013060193
#> 2               0.02618614                 0.0102981933
#> 3              -0.10998772                 0.0096872090
#> 4              -0.29960088                -0.0007205544
#> 5              -0.33078317                 0.0189097688
#> 6               0.89746122                 0.0110337423
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                       0.00000000
#> 2                 0.6666667                       0.84444444
#> 3                 0.5333333                       0.48888889
#> 4                 0.0000000                       0.04444444
#> 5                 0.0000000                       0.00000000
#> 6                 0.0000000                       0.00000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.4000000
#> 3                            0.5777778
#> 4                            0.6222222
#> 5                            0.4000000
#> 6                            0.0000000
```

We can now plot this results as trait networks using the
**plotNetworks**. Let’s first plot the phylogenetic part of the variance
and covariance.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  1                  2                  1                  2

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.39 
    #> 2     ---     4   -0.33 
    #> 3     ---     4   0.9

Let’s now plot the non-phylogenetic variances and covariances.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  1                  1                  0                  0

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.5

We observe how much of the variance and covariance in the phylo_G1
traits is related to the phylogeny, and that the variances and
covariances in nonPhylo_G1 traits are not related to the phylogeny, as
expected.

## Compute variances and covariances for simulated data including an enviornmental variable

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
    #> 1 nonPhylo_G1_trait1                 100             0.2714437
    #> 2 nonPhylo_G1_trait2                 100             0.2258600
    #> 3    phylo_G1_trait1                 100             0.9665661
    #> 4    phylo_G1_trait2                 100             0.9515041
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1                0.72855631                             0
    #> 2                0.75724742                             0
    #> 3                0.01701295                             0
    #> 4                0.01751488                             0
    #>   environmental_variable non_attributed_phylogenetic_variance
    #> 1           phylo_G1_env                            0.2374249
    #> 2           phylo_G1_env                            0.1979330
    #> 3           phylo_G1_env                            0.3771059
    #> 4           phylo_G1_env                            0.4743021
    #>   environmental_phylogenetic_variance labile_environmental_variance
    #> 1                          0.03401876                    0.03482793
    #> 2                          0.02792699                    0.02875520
    #> 3                          0.58946020                    0.01633438
    #> 4                          0.47720200                    0.03387704
    #>   residual_variance p_value_non_attributed_phylogenetic_variance
    #> 1        0.69372838                                            0
    #> 2        0.74538483                                            0
    #> 3        0.01709951                                            0
    #> 4        0.01461886                                            0
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
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.83824250
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.26135716
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.13659127
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.13000302
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.04270588
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90304397
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.19549236                  0.642750139
    #> 2              -0.25883810                 -0.002519060
    #> 3              -0.17593092                  0.039339643
    #> 4              -0.12442893                 -0.005574089
    #> 5              -0.05509789                  0.012392017
    #> 6               0.88239316                  0.020650815
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                              0.0
    #> 2                 0.0000000                              0.0
    #> 3                 0.0000000                              0.0
    #> 4                 0.4888889                              0.4
    #> 5                 0.7555556                              0.4
    #> 6                 0.0000000                              0.0
    #>   p_value_non_phylogenetic_correlation environmental_variable
    #> 1                           0.00000000           phylo_G1_env
    #> 2                           0.88888889           phylo_G1_env
    #> 3                           0.08888889           phylo_G1_env
    #> 4                           0.84444444           phylo_G1_env
    #> 5                           0.57777778           phylo_G1_env
    #> 6                           0.00000000           phylo_G1_env
    #>   non_attributed_phylogenetic_correlation
    #> 1                            1.480543e-01
    #> 2                           -1.276528e-01
    #> 3                           -1.434944e-01
    #> 4                           -2.479039e-05
    #> 5                           -5.288616e-02
    #> 6                            3.353705e-01
    #>   environmental_phylogenetic_correlation labile_environmental_correlation
    #> 1                            0.047438047                      0.026697653
    #> 2                           -0.131185269                      0.006251765
    #> 3                           -0.032436481                      0.020559955
    #> 4                           -0.124404144                      0.001914865
    #> 5                           -0.002211731                      0.004907669
    #> 6                            0.547022637                      0.014716695
    #>   residual_correlation p_value_non_attributed_phylogenetic_correlation
    #> 1          0.616052486                                      0.00000000
    #> 2         -0.008770826                                      0.08888889
    #> 3          0.018779688                                      0.00000000
    #> 4         -0.007488954                                      0.97777778
    #> 5          0.007484349                                      0.31111111
    #> 6          0.005934120                                      0.00000000
    #>   p_value_environmental_phylogenetic_correlation
    #> 1                                     0.00000000
    #> 2                                     0.04444444
    #> 3                                     0.66666667
    #> 4                                     0.00000000
    #> 5                                     0.93333333
    #> 6                                     0.00000000
    #>   p_value_labile_environmental_correlation p_value_residual_correlation
    #> 1                                0.1777778                    0.0000000
    #> 2                                0.7555556                    0.5777778
    #> 3                                0.2666667                    0.3111111
    #> 4                                0.7111111                    0.7111111
    #> 5                                0.6666667                    0.7111111
    #> 6                                0.0000000                    0.0000000

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2714437
#> 2 nonPhylo_G1_trait2                 100             0.2258600
#> 3    phylo_G1_trait1                 100             0.9665661
#> 4    phylo_G1_trait2                 100             0.9515041
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1                0.72855631                             0
#> 2                0.75724742                             0
#> 3                0.01701295                             0
#> 4                0.01751488                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.2374249
#> 2           phylo_G1_env                            0.1979330
#> 3           phylo_G1_env                            0.3771059
#> 4           phylo_G1_env                            0.4743021
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.03401876                    0.03482793
#> 2                          0.02792699                    0.02875520
#> 3                          0.58946020                    0.01633438
#> 4                          0.47720200                    0.03387704
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.69372838                                            0
#> 2        0.74538483                                            0
#> 3        0.01709951                                            0
#> 4        0.01461886                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.83824250
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.26135716
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.13659127
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.13000302
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.04270588
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90304397
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.19549236                  0.642750139
#> 2              -0.25883810                 -0.002519060
#> 3              -0.17593092                  0.039339643
#> 4              -0.12442893                 -0.005574089
#> 5              -0.05509789                  0.012392017
#> 6               0.88239316                  0.020650815
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                              0.0
#> 2                 0.0000000                              0.0
#> 3                 0.0000000                              0.0
#> 4                 0.4888889                              0.4
#> 5                 0.7555556                              0.4
#> 6                 0.0000000                              0.0
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                           0.00000000           phylo_G1_env
#> 2                           0.88888889           phylo_G1_env
#> 3                           0.08888889           phylo_G1_env
#> 4                           0.84444444           phylo_G1_env
#> 5                           0.57777778           phylo_G1_env
#> 6                           0.00000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                            1.480543e-01
#> 2                           -1.276528e-01
#> 3                           -1.434944e-01
#> 4                           -2.479039e-05
#> 5                           -5.288616e-02
#> 6                            3.353705e-01
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.047438047                      0.026697653
#> 2                           -0.131185269                      0.006251765
#> 3                           -0.032436481                      0.020559955
#> 4                           -0.124404144                      0.001914865
#> 5                           -0.002211731                      0.004907669
#> 6                            0.547022637                      0.014716695
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.616052486                                      0.00000000
#> 2         -0.008770826                                      0.08888889
#> 3          0.018779688                                      0.00000000
#> 4         -0.007488954                                      0.97777778
#> 5          0.007484349                                      0.31111111
#> 6          0.005934120                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                     0.00000000
#> 2                                     0.04444444
#> 3                                     0.66666667
#> 4                                     0.00000000
#> 5                                     0.93333333
#> 6                                     0.00000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                                0.1777778                    0.0000000
#> 2                                0.7555556                    0.5777778
#> 3                                0.2666667                    0.3111111
#> 4                                0.7111111                    0.7111111
#> 5                                0.6666667                    0.7111111
#> 6                                0.0000000                    0.0000000
```

Let’s plot the phylogenetic environmental variances and covariances as
an example.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  0                  0                  1                  1

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.55

We can see how part of the pattern of variances and covariances of the
traits phylo_G1 are attributed to the phylogenetically conserved effect
of the environmental variable phylo_G1_env.

## Impute missing values

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
    #> 1    phylo_G1_env                  80             0.9722607
    #> 2 phylo_G1_trait1                  80             0.9773423
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1                0.02773932                             0
    #> 2                0.02360152                             0
    #> [1] "Covariance results:"
    #>        trait_1         trait_2 number_observations total_correlation
    #> 1 phylo_G1_env phylo_G1_trait1                  80         0.7581911
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1                0.7347814                   0.02340977
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
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7124264
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7076517
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7109043
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.6981408
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.710597
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7137866
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.726343
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.72128
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7247107
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7197685
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7023334
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.721973
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7003642
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.696556
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7376661
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7133629
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7244556
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7313482
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7209711
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7285549
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7305749
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.708448
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7290773
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7281112
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7152948
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.716226
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7070215
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7114473
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.708177
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA     NRMSE
    #> 1 phylo_G1_trait1 100    80   20 0.7192802
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1        1.843047
#> 2    t10       -1.494269
#> 3   t100       -1.682892
#> 4    t11        3.447466
#> 5    t12       -2.662426
#> 6    t13       -4.460927
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
#> -6.4114 -2.1948 -0.7595 -0.7744  0.6813  3.5858      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -6.4114 -1.8416 -0.7951 -0.7800  0.5613  3.5858
```
