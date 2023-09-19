
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
#> 1     t1      -3.7569593       -6.293649   -2.1433150       3.3470698
#> 2    t10      -5.3021602       -5.596330   -6.2135112       0.9596298
#> 3   t100      -1.4687763       -1.999095   -3.5745067       1.1963119
#> 4    t11      -0.8281509       -2.048837   -2.8759795       0.4069202
#> 5    t12      -1.7592195       -2.399128   -1.0967148      -1.0517331
#> 6    t13      -2.3414466       -3.606760    0.9329666      -6.4318589
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1       2.5945102    1.3114326          0.5936224         0.72494965
#> 2       1.3549047   -0.2784216          0.5752437         0.67062138
#> 3       2.6356911   -0.4850631         -0.5991368         0.02737421
#> 4       1.7768562   -2.2877624          0.5190405        -0.44824448
#> 5      -0.4930878   -0.5419400          0.2759576        -0.08632709
#> 6      -6.8932125   -5.4070417          0.3729064         0.79230235
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1      2.09628461         -0.4815900          0.1181552      -0.1475090
#> 2     -0.08074836          0.1747601          0.4718272       0.8832021
#> 3     -0.48096823          0.7449276          0.3313209       0.8577608
#> 4     -0.23351860          1.4039725          1.6839296       1.1617792
#> 5      0.07897095          0.0660290         -0.2511019       0.2375623
#> 6      0.34733121          0.6941880          1.5484659       1.0942025
```

We can use the function *plotData* to plot some simulated traits. Let’s
plot the phylo_G1 traits and environmental variable.

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

In this plot we can se how the simulated variances and covariances are
phylogenetically conserved. Let’s now plot the nonPhylo_G1 traits and
environmental variable.

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
    #> 1 nonPhylo_G1_trait1                 100             0.3908425
    #> 2 nonPhylo_G1_trait2                 100             0.2822254
    #> 3    phylo_G1_trait1                 100             0.9942257
    #> 4    phylo_G1_trait2                 100             0.9778130
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.609157515                             0
    #> 2               0.692214072                             0
    #> 3               0.008632071                             0
    #> 4               0.010285776                             0
    #> [1] "Covariance results:"
    #>              trait_1            trait_2 number_observations total_correlation
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100         0.8092536
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        -0.3009209
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        -0.1233403
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        -0.1409152
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        -0.1877545
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100         0.8940046
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1                0.3508199                  0.458433687
    #> 2               -0.3392546                  0.038333713
    #> 3               -0.1553548                  0.032014522
    #> 4               -0.1937654                  0.052850265
    #> 5               -0.1940912                  0.006336718
    #> 6                0.8894988                  0.004505833
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                       0.00000000
    #> 2                 0.0000000                       0.00000000
    #> 3                 0.3111111                       0.26666667
    #> 4                 0.1777778                       0.13333333
    #> 5                 0.0000000                       0.04444444
    #> 6                 0.0000000                       0.00000000
    #>   p_value_non_phylogenetic_correlation
    #> 1                           0.00000000
    #> 2                           0.00000000
    #> 3                           0.31111111
    #> 4                           0.08888889
    #> 5                           0.75555556
    #> 6                           0.08888889

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.3908425
#> 2 nonPhylo_G1_trait2                 100             0.2822254
#> 3    phylo_G1_trait1                 100             0.9942257
#> 4    phylo_G1_trait2                 100             0.9778130
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.609157515                             0
#> 2               0.692214072                             0
#> 3               0.008632071                             0
#> 4               0.010285776                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100         0.8092536
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        -0.3009209
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        -0.1233403
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100        -0.1409152
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        -0.1877545
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100         0.8940046
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1                0.3508199                  0.458433687
#> 2               -0.3392546                  0.038333713
#> 3               -0.1553548                  0.032014522
#> 4               -0.1937654                  0.052850265
#> 5               -0.1940912                  0.006336718
#> 6                0.8894988                  0.004505833
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                       0.00000000
#> 2                 0.0000000                       0.00000000
#> 3                 0.3111111                       0.26666667
#> 4                 0.1777778                       0.13333333
#> 5                 0.0000000                       0.04444444
#> 6                 0.0000000                       0.00000000
#>   p_value_non_phylogenetic_correlation
#> 1                           0.00000000
#> 2                           0.00000000
#> 3                           0.31111111
#> 4                           0.08888889
#> 5                           0.75555556
#> 6                           0.08888889
```

We can now plot this results as trait networks using the
**plotNetworks**. Let’s first plot the phylogenetic part of the variance
and covariance.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  2                  1                  2                  1

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.35 
    #> 1     ---     3   -0.34 
    #> 3     ---     4   0.89

Let’s now plot the non-phylogenetic variances and covariances.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  1                  1                  0                  0

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

    #> From     To  Weight
    #> 1     ---     2   0.46

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
    #> 1 nonPhylo_G1_trait1                 100             0.2641821
    #> 2 nonPhylo_G1_trait2                 100             0.2746948
    #> 3    phylo_G1_trait1                 100             0.9784213
    #> 4    phylo_G1_trait2                 100             0.9803976
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.735817928                             0
    #> 2               0.681320627                             0
    #> 3               0.007384811                             0
    #> 4              -0.005388570                             0
    #>   environmental_variable non_attributed_phylogenetic_variance
    #> 1           phylo_G1_env                            0.2560755
    #> 2           phylo_G1_env                            0.2519098
    #> 3           phylo_G1_env                            0.3097768
    #> 4           phylo_G1_env                            0.3761429
    #>   environmental_phylogenetic_variance labile_environmental_variance
    #> 1                         0.008106598                   0.030402795
    #> 2                         0.022784981                   0.074638845
    #> 3                         0.668644476                   0.006238181
    #> 4                         0.604254765                   0.005711434
    #>   residual_variance p_value_non_attributed_phylogenetic_variance
    #> 1        0.70541513                                            0
    #> 2        0.65066638                                            0
    #> 3        0.01534052                                            0
    #> 4        0.01389095                                            0
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
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.72466452
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.14735752
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.08126023
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.09062170
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.08491360
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90386725
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.13033475                   0.59432977
    #> 2              -0.13449889                  -0.01285863
    #> 3              -0.06461326                  -0.01664697
    #> 4              -0.12828923                   0.03766754
    #> 5              -0.12797861                   0.04306501
    #> 6               0.88984490                   0.01402235
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                        0.3555556
    #> 2                 0.1777778                        0.2222222
    #> 3                 0.2666667                        0.3555556
    #> 4                 0.2666667                        0.1777778
    #> 5                 0.6222222                        0.4000000
    #> 6                 0.0000000                        0.0000000
    #>   p_value_non_phylogenetic_correlation environmental_variable
    #> 1                           0.00000000           phylo_G1_env
    #> 2                           0.62222222           phylo_G1_env
    #> 3                           0.44444444           phylo_G1_env
    #> 4                           0.04444444           phylo_G1_env
    #> 5                           0.08888889           phylo_G1_env
    #> 6                           0.00000000           phylo_G1_env
    #>   non_attributed_phylogenetic_correlation
    #> 1                              0.12576886
    #> 2                             -0.05384316
    #> 3                             -0.01124836
    #> 4                             -0.12134618
    #> 5                             -0.03592373
    #> 6                              0.24954112
    #>   environmental_phylogenetic_correlation labile_environmental_correlation
    #> 1                            0.004565896                      0.033535549
    #> 2                           -0.080655725                     -0.005766677
    #> 3                           -0.053364900                     -0.011276959
    #> 4                           -0.006943047                      0.010304086
    #> 5                           -0.092054880                      0.021332496
    #> 6                            0.640303778                      0.004849270
    #>   residual_correlation p_value_non_attributed_phylogenetic_correlation
    #> 1          0.560794219                                       0.3111111
    #> 2         -0.007091958                                       0.1777778
    #> 3         -0.005370012                                       0.8888889
    #> 4          0.027363450                                       0.0000000
    #> 5          0.021732509                                       0.7111111
    #> 6          0.009173075                                       0.0000000
    #>   p_value_environmental_phylogenetic_correlation
    #> 1                                      0.8444444
    #> 2                                      0.4888889
    #> 3                                      0.5333333
    #> 4                                      0.9777778
    #> 5                                      0.4888889
    #> 6                                      0.0000000
    #>   p_value_labile_environmental_correlation p_value_residual_correlation
    #> 1                                0.4000000                   0.00000000
    #> 2                                0.6666667                   0.57777778
    #> 3                                0.1777778                   0.80000000
    #> 4                                0.3555556                   0.04444444
    #> 5                                0.1777778                   0.35555556
    #> 6                                0.2222222                   0.00000000

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2641821
#> 2 nonPhylo_G1_trait2                 100             0.2746948
#> 3    phylo_G1_trait1                 100             0.9784213
#> 4    phylo_G1_trait2                 100             0.9803976
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.735817928                             0
#> 2               0.681320627                             0
#> 3               0.007384811                             0
#> 4              -0.005388570                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.2560755
#> 2           phylo_G1_env                            0.2519098
#> 3           phylo_G1_env                            0.3097768
#> 4           phylo_G1_env                            0.3761429
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                         0.008106598                   0.030402795
#> 2                         0.022784981                   0.074638845
#> 3                         0.668644476                   0.006238181
#> 4                         0.604254765                   0.005711434
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.70541513                                            0
#> 2        0.65066638                                            0
#> 3        0.01534052                                            0
#> 4        0.01389095                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.72466452
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.14735752
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.08126023
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.09062170
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.08491360
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90386725
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.13033475                   0.59432977
#> 2              -0.13449889                  -0.01285863
#> 3              -0.06461326                  -0.01664697
#> 4              -0.12828923                   0.03766754
#> 5              -0.12797861                   0.04306501
#> 6               0.88984490                   0.01402235
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                        0.3555556
#> 2                 0.1777778                        0.2222222
#> 3                 0.2666667                        0.3555556
#> 4                 0.2666667                        0.1777778
#> 5                 0.6222222                        0.4000000
#> 6                 0.0000000                        0.0000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                           0.00000000           phylo_G1_env
#> 2                           0.62222222           phylo_G1_env
#> 3                           0.44444444           phylo_G1_env
#> 4                           0.04444444           phylo_G1_env
#> 5                           0.08888889           phylo_G1_env
#> 6                           0.00000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.12576886
#> 2                             -0.05384316
#> 3                             -0.01124836
#> 4                             -0.12134618
#> 5                             -0.03592373
#> 6                              0.24954112
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.004565896                      0.033535549
#> 2                           -0.080655725                     -0.005766677
#> 3                           -0.053364900                     -0.011276959
#> 4                           -0.006943047                      0.010304086
#> 5                           -0.092054880                      0.021332496
#> 6                            0.640303778                      0.004849270
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.560794219                                       0.3111111
#> 2         -0.007091958                                       0.1777778
#> 3         -0.005370012                                       0.8888889
#> 4          0.027363450                                       0.0000000
#> 5          0.021732509                                       0.7111111
#> 6          0.009173075                                       0.0000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.8444444
#> 2                                      0.4888889
#> 3                                      0.5333333
#> 4                                      0.9777778
#> 5                                      0.4888889
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                                0.4000000                   0.00000000
#> 2                                0.6666667                   0.57777778
#> 3                                0.1777778                   0.80000000
#> 4                                0.3555556                   0.04444444
#> 5                                0.1777778                   0.35555556
#> 6                                0.2222222                   0.00000000
```

Let’s plot the phylogenetic environmental variances and covariances as
an example.

    #> [1] "Node degree:"
    #> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
    #>                  0                  0                  1                  1

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.64

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
    #> 1    phylo_G1_env                  80             0.9809437
    #> 2 phylo_G1_trait1                  80             0.9849273
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1                0.01905628                             0
    #> 2                0.01458124                             0
    #> [1] "Covariance results:"
    #>        trait_1         trait_2 number_observations total_correlation
    #> 1 phylo_G1_env phylo_G1_trait1                  80         0.7816955
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1                0.7665198                   0.01517574
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                         0                                0
    #>   p_value_non_phylogenetic_correlation
    #> 1                                    0
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.907146
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.913054
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.915361
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.936757
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.923334
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.904445
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.93773
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.940003
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.961822
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.987204
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.964711
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.974657
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.954199
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.929149
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.906689
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.946855
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.982741
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.92683
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.997341
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.929291
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.977805
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.943079
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.959473
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.927283
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.98815
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.996646
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.951155
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.924244
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.975207
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.976516
    #>                                                                                                                                                                                              Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, phylo_G1_env

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1      -3.7569593
#> 2    t10      -5.3021602
#> 3   t100      -1.4687763
#> 4    t11      -0.8281509
#> 5    t12      -1.7592195
#> 6    t13      -2.3414466
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
#> -9.3806 -4.3173 -2.8940 -1.9213  0.7535  5.1229      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -9.381  -4.137  -2.544  -1.951  -0.538   5.123
```
