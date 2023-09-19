
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
#> 1     t1      2.83569287       3.5784519   5.15925627      -0.9055044
#> 2    t10      0.01988534       0.3293684   1.78411002      -3.6064855
#> 3   t100      0.30868880       1.1886256  -1.94854582      -4.2350166
#> 4    t11      1.38680507      -0.3226524   0.04940396      -2.3785741
#> 5    t12      0.47420341       1.4549272  -1.60981599      -3.9457055
#> 6    t13      0.87811755      -0.2289799   1.82860671       1.9071952
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1       0.3714268   -1.8444343         -0.4471672         -0.7720925
#> 2      -1.1111146   -4.2793764          1.3149870          1.2467222
#> 3      -2.3824668   -3.7698756         -0.7394864          0.0991954
#> 4      -1.3602870   -3.2604774          1.3968993          1.2233685
#> 5      -2.5368309   -3.3527166          0.2102929          0.6311535
#> 6       1.7904643   -0.7320168         -2.0835196         -2.1369207
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1     -1.64095407         -0.6458856        -0.65334877     -0.82827619
#> 2      1.42099324         -0.6283133        -0.04520841     -0.60510485
#> 3     -0.22839588          0.5680749         0.51480651     -0.01322936
#> 4      0.63844099          0.6960582         0.24112295      0.85966124
#> 5     -0.02877568         -0.3140729         0.07506147      0.08764462
#> 6     -1.28750430          0.4565116         0.46011956     -0.99582493
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
    #> 1 nonPhylo_G1_trait1                 100             0.1216253
    #> 2 nonPhylo_G1_trait2                 100             0.2444461
    #> 3    phylo_G1_trait1                 100             0.9776827
    #> 4    phylo_G1_trait2                 100             0.9604399
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.878374674                             0
    #> 2               0.751036373                             0
    #> 3               0.019350253                             0
    #> 4              -0.008784689                             0
    #> [1] "Covariance results:"
    #>              trait_1            trait_2 number_observations total_correlation
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100       0.854371267
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       0.022845219
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       0.337260502
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100      -0.028439522
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       0.001228066
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100       0.863613407
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.08103834                  0.773332923
    #> 2               0.02039136                  0.002453861
    #> 3               0.35428683                 -0.017026323
    #> 4              -0.05280336                  0.024363843
    #> 5               0.01350667                 -0.012278607
    #> 6               0.82155443                  0.042058973
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                        0.3111111
    #> 2                 0.8000000                        0.8444444
    #> 3                 0.0000000                        0.0000000
    #> 4                 0.8444444                        0.6666667
    #> 5                 0.8444444                        0.8888889
    #> 6                 0.0000000                        0.0000000
    #>   p_value_non_phylogenetic_correlation
    #> 1                            0.0000000
    #> 2                            0.7555556
    #> 3                            0.3111111
    #> 4                            0.1777778
    #> 5                            0.2666667
    #> 6                            0.0000000

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.1216253
#> 2 nonPhylo_G1_trait2                 100             0.2444461
#> 3    phylo_G1_trait1                 100             0.9776827
#> 4    phylo_G1_trait2                 100             0.9604399
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.878374674                             0
#> 2               0.751036373                             0
#> 3               0.019350253                             0
#> 4              -0.008784689                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100       0.854371267
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       0.022845219
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       0.337260502
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100      -0.028439522
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       0.001228066
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100       0.863613407
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.08103834                  0.773332923
#> 2               0.02039136                  0.002453861
#> 3               0.35428683                 -0.017026323
#> 4              -0.05280336                  0.024363843
#> 5               0.01350667                 -0.012278607
#> 6               0.82155443                  0.042058973
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                        0.3111111
#> 2                 0.8000000                        0.8444444
#> 3                 0.0000000                        0.0000000
#> 4                 0.8444444                        0.6666667
#> 5                 0.8444444                        0.8888889
#> 6                 0.0000000                        0.0000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.7555556
#> 3                            0.3111111
#> 4                            0.1777778
#> 5                            0.2666667
#> 6                            0.0000000
```

We can now plot this results as trait networks using the
**plotNetworks**. Let’s first plot the phylogenetic part of the variance
and covariance.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  1                  0                  1                  2

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     4   0.35 
    #> 3     ---     4   0.82

Let’s now plot the non-phylogenetic variances and covariances.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  1                  1                  0                  0

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.77

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
    #> 1 nonPhylo_G1_trait1                 100            0.09822224
    #> 2 nonPhylo_G1_trait2                 100            0.23110807
    #> 3    phylo_G1_trait1                 100            0.96849176
    #> 4    phylo_G1_trait2                 100            0.93944726
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.901777757                             0
    #> 2               0.762852973                             0
    #> 3               0.007001850                             0
    #> 4               0.005015341                             0
    #>   environmental_variable non_attributed_phylogenetic_variance
    #> 1           phylo_G1_env                            0.0955155
    #> 2           phylo_G1_env                            0.2056510
    #> 3           phylo_G1_env                            0.4174514
    #> 4           phylo_G1_env                            0.4011485
    #>   environmental_phylogenetic_variance labile_environmental_variance
    #> 1                         0.002706738                    0.03750646
    #> 2                         0.025457029                    0.06767015
    #> 3                         0.551040397                    0.01348227
    #> 4                         0.538298792                    0.02474695
    #>   residual_variance p_value_non_attributed_phylogenetic_variance
    #> 1        0.86427130                                            0
    #> 2        0.70122178                                            0
    #> 3        0.01802598                                            0
    #> 4        0.03580579                                            0
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
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.85455142
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.22119861
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.08292022
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.13356124
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        0.02871787
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.89132402
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.06389783                  0.790653586
    #> 2               0.25175172                 -0.030553111
    #> 3               0.07298193                  0.009938282
    #> 4              -0.14277353                  0.009212297
    #> 5               0.02045923                  0.008258640
    #> 6               0.86789993                  0.023424092
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                0.00000000                       0.00000000
    #> 2                0.04444444                       0.00000000
    #> 3                0.35555556                       0.31111111
    #> 4                0.04444444                       0.08888889
    #> 5                0.84444444                       0.97777778
    #> 6                0.00000000                       0.00000000
    #>   p_value_non_phylogenetic_correlation environmental_variable
    #> 1                            0.0000000           phylo_G1_env
    #> 2                            0.4000000           phylo_G1_env
    #> 3                            0.8444444           phylo_G1_env
    #> 4                            0.7555556           phylo_G1_env
    #> 5                            0.9777778           phylo_G1_env
    #> 6                            0.0000000           phylo_G1_env
    #>   non_attributed_phylogenetic_correlation
    #> 1                              0.06044406
    #> 2                              0.02493880
    #> 3                             -0.01485658
    #> 4                             -0.06908180
    #> 5                             -0.08376587
    #> 6                              0.27926576
    #>   environmental_phylogenetic_correlation labile_environmental_correlation
    #> 1                             0.00345377                      0.037877442
    #> 2                             0.22681292                     -0.031456348
    #> 3                             0.08783852                     -0.012218320
    #> 4                            -0.07369173                     -0.001754469
    #> 5                             0.10422510                     -0.012466113
    #> 6                             0.58863417                      0.012702531
    #>   residual_correlation p_value_non_attributed_phylogenetic_correlation
    #> 1         0.7527761447                                       0.0000000
    #> 2         0.0009032369                                       0.6222222
    #> 3         0.0221566018                                       0.6222222
    #> 4         0.0109667660                                       0.1777778
    #> 5         0.0207247534                                       0.2666667
    #> 6         0.0107215613                                       0.0000000
    #>   p_value_environmental_phylogenetic_correlation
    #> 1                                      0.6222222
    #> 2                                      0.0000000
    #> 3                                      0.4000000
    #> 4                                      0.5333333
    #> 5                                      0.1777778
    #> 6                                      0.0000000
    #>   p_value_labile_environmental_correlation p_value_residual_correlation
    #> 1                                0.2666667                   0.00000000
    #> 2                                0.2222222                   0.93333333
    #> 3                                0.4000000                   0.53333333
    #> 4                                0.9777778                   0.48888889
    #> 5                                0.4888889                   0.53333333
    #> 6                                0.8000000                   0.08888889

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100            0.09822224
#> 2 nonPhylo_G1_trait2                 100            0.23110807
#> 3    phylo_G1_trait1                 100            0.96849176
#> 4    phylo_G1_trait2                 100            0.93944726
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.901777757                             0
#> 2               0.762852973                             0
#> 3               0.007001850                             0
#> 4               0.005015341                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.0955155
#> 2           phylo_G1_env                            0.2056510
#> 3           phylo_G1_env                            0.4174514
#> 4           phylo_G1_env                            0.4011485
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                         0.002706738                    0.03750646
#> 2                         0.025457029                    0.06767015
#> 3                         0.551040397                    0.01348227
#> 4                         0.538298792                    0.02474695
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.86427130                                            0
#> 2        0.70122178                                            0
#> 3        0.01802598                                            0
#> 4        0.03580579                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.85455142
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.22119861
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.08292022
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.13356124
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        0.02871787
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.89132402
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.06389783                  0.790653586
#> 2               0.25175172                 -0.030553111
#> 3               0.07298193                  0.009938282
#> 4              -0.14277353                  0.009212297
#> 5               0.02045923                  0.008258640
#> 6               0.86789993                  0.023424092
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.00000000
#> 2                0.04444444                       0.00000000
#> 3                0.35555556                       0.31111111
#> 4                0.04444444                       0.08888889
#> 5                0.84444444                       0.97777778
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                            0.0000000           phylo_G1_env
#> 2                            0.4000000           phylo_G1_env
#> 3                            0.8444444           phylo_G1_env
#> 4                            0.7555556           phylo_G1_env
#> 5                            0.9777778           phylo_G1_env
#> 6                            0.0000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.06044406
#> 2                              0.02493880
#> 3                             -0.01485658
#> 4                             -0.06908180
#> 5                             -0.08376587
#> 6                              0.27926576
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                             0.00345377                      0.037877442
#> 2                             0.22681292                     -0.031456348
#> 3                             0.08783852                     -0.012218320
#> 4                            -0.07369173                     -0.001754469
#> 5                             0.10422510                     -0.012466113
#> 6                             0.58863417                      0.012702531
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1         0.7527761447                                       0.0000000
#> 2         0.0009032369                                       0.6222222
#> 3         0.0221566018                                       0.6222222
#> 4         0.0109667660                                       0.1777778
#> 5         0.0207247534                                       0.2666667
#> 6         0.0107215613                                       0.0000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.6222222
#> 2                                      0.0000000
#> 3                                      0.4000000
#> 4                                      0.5333333
#> 5                                      0.1777778
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                                0.2666667                   0.00000000
#> 2                                0.2222222                   0.93333333
#> 3                                0.4000000                   0.53333333
#> 4                                0.9777778                   0.48888889
#> 5                                0.4888889                   0.53333333
#> 6                                0.8000000                   0.08888889
```

Let’s plot the phylogenetic environmental variances and covariances as
an example.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  0                  0                  1                  1

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.59

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
    #> 1    phylo_G1_env                  80             0.9794678
    #> 2 phylo_G1_trait1                  80             0.9616898
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1                0.02053223                             0
    #> 2                0.02246987                             0
    #> [1] "Covariance results:"
    #>        trait_1         trait_2 number_observations total_correlation
    #> 1 phylo_G1_env phylo_G1_trait1                  80         0.6276755
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1                0.6054601                   0.02221536
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
    #> 1 phylo_G1_trait1 100    80   20 1.223724
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.223307
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.238298
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.211371
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.19276
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.209378
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.247754
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.237416
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.223641
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.230015
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.246626
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.220412
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.209345
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.203584
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.205562
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.202263
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.208889
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.224326
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.208624
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.245801
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.205815
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.234414
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.228487
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA NRMSE
    #> 1 phylo_G1_trait1 100    80   20  1.23
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.209134
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.215634
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.212565
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.189259
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.255474
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.22282
    #>                                                                                                                                                                                                             Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1      2.83569287
#> 2    t10      0.01988534
#> 3   t100      0.30868880
#> 4    t11      2.61922370
#> 5    t12      0.47420341
#> 6    t13      0.87811755
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
#> -6.4282 -0.9339  0.7182  0.3193  2.0054  4.8047      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -6.4282 -1.2303  0.6649  0.2981  2.0314  4.8047
```
