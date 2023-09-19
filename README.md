
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
#> 1     t1       3.5950643      4.01968414    4.8263179       -1.067000
#> 2    t10       1.3589125     -0.02661093   -0.2651821        6.666912
#> 3   t100       5.3216593      3.74118087    0.5330322        4.489373
#> 4    t11       0.1944281     -3.11042498   -2.3143392        1.876429
#> 5    t12      -3.2008335     -4.78106964   -5.5359417        1.965355
#> 6    t13      -0.1135998     -4.36345714   -5.7605852       -8.246929
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1      -2.7012115     3.836080         -0.2060912          0.2033346
#> 2       6.3781491     5.226321          0.1922425         -0.6550333
#> 3       3.6170939     3.797612          0.3965412         -0.4525513
#> 4      -1.0033583     1.071417         -1.6928618         -2.3302311
#> 5      -0.7424521     1.119069         -0.8579106         -1.0073838
#> 6     -10.0518107    -9.209030         -0.7328977         -1.2680203
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1       0.3264828         -0.8245313        -0.31402669      0.05733168
#> 2       0.1370826          0.6436855         0.28169871      1.27852633
#> 3      -0.6863814          0.6385627        -0.36073214      0.02513247
#> 4      -1.8481146          0.6304362         0.97443414     -0.16365350
#> 5      -0.6114646         -0.2389539         0.02798378     -0.23614438
#> 6      -2.0487201         -0.2112777        -0.96390854     -1.14585038
```

We can use the function *plotData* to plot some simulated traits. Let’s
plot the phylo_G1 traits and environmental variable.

``` r
plotData(variables = c("phylo_G1_trait1", "phylo_G1_trait2", "phylo_G1_env"),
         dataset = simulated.data$data,
         phylogeny = simulated.data$phylogeny,
         terminal.taxon = "animal" 
          )
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

In this plot we can se how the simulated variances and covariances are
phylogenetically conserved. Let’s now plot the nonPhylo_G1 traits and
environmental variable.

``` r
plotData(variables = c("nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "nonPhylo_G1_env"),
         dataset = simulated.data$data,
         phylogeny = simulated.data$phylogeny,
         terminal.taxon = "animal" 
          )
```

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

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.3101489
#> 2 nonPhylo_G1_trait2                 100             0.2407941
#> 3    phylo_G1_trait1                 100             0.9787773
#> 4    phylo_G1_trait2                 100             0.9834053
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.689851144                             0
#> 2               0.733807211                             0
#> 3               0.011635448                             0
#> 4               0.006943518                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.84882897
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.11021249
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.03122048
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.27215798
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        0.13374918
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.88897133
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.16314171                 0.6856872543
#> 2               0.10812887                 0.0020836162
#> 3              -0.05437862                 0.0231581398
#> 4              -0.28143250                 0.0092745269
#> 5               0.13284147                 0.0009077053
#> 6               0.86809799                 0.0208733391
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                        0.0000000
#> 2                 0.5333333                        0.5333333
#> 3                 0.5777778                        0.4000000
#> 4                 0.0000000                        0.0000000
#> 5                 0.1333333                        0.1777778
#> 6                 0.0000000                        0.0000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.9333333
#> 3                            0.0000000
#> 4                            0.4888889
#> 5                            0.9333333
#> 6                            0.0000000
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
            correlation_threshold = 0.3,
            label_size = 1)
#> [1] "Node degree:"
#> nonPhylo_G1_trait1 nonPhylo_G1_trait2    phylo_G1_trait1    phylo_G1_trait2 
#>                  0                  0                  1                  1
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

    #> From     To  Weight
    #> 3     ---     4   0.87

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
    #> 1     ---     2   0.69

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
#> 1 nonPhylo_G1_trait1                 100             0.1940665
#> 2 nonPhylo_G1_trait2                 100             0.3108555
#> 3    phylo_G1_trait1                 100             0.9748093
#> 4    phylo_G1_trait2                 100             0.9842086
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.805933481                             0
#> 2               0.670072881                             0
#> 3               0.009314317                             0
#> 4               0.004464561                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1695252
#> 2           phylo_G1_env                            0.2787971
#> 3           phylo_G1_env                            0.3838660
#> 4           phylo_G1_env                            0.3260598
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.02454136                   0.025041681
#> 2                          0.03205846                   0.035496272
#> 3                          0.59094328                   0.015734893
#> 4                          0.65814885                   0.005949253
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1       0.780891800                                            0
#> 2       0.653648193                                            0
#> 3       0.009455851                                            0
#> 4       0.009842140                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100       0.761780431
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100      -0.071022997
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100      -0.005413472
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       0.271162837
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       0.076093279
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100       0.923916821
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.05372657                  0.708053859
#> 2              -0.07466198                  0.003638987
#> 3              -0.01872596                  0.013312492
#> 4               0.28854708                 -0.017384243
#> 5               0.07291982                  0.003173455
#> 6               0.91376483                  0.010151992
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.22222222
#> 2                0.71111111                       0.57777778
#> 3                0.97777778                       0.84444444
#> 4                0.04444444                       0.04444444
#> 5                0.53333333                       0.53333333
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                            0.0000000           phylo_G1_env
#> 2                            0.6666667           phylo_G1_env
#> 3                            0.5333333           phylo_G1_env
#> 4                            0.5777778           phylo_G1_env
#> 5                            0.8888889           phylo_G1_env
#> 6                            0.0000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.04914878
#> 2                              0.01935421
#> 3                              0.06463409
#> 4                              0.12048005
#> 5                              0.10480187
#> 6                              0.27282310
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.004577788                     0.0240248096
#> 2                           -0.094016191                    -0.0005461307
#> 3                           -0.083360050                     0.0064731785
#> 4                            0.168067029                    -0.0180825270
#> 5                           -0.031882041                     0.0006665353
#> 6                            0.640941731                     0.0067273769
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1         0.6840290489                                       0.2222222
#> 2         0.0041851178                                       0.7111111
#> 3         0.0068393133                                       0.2666667
#> 4         0.0006982836                                       0.0000000
#> 5         0.0025069194                                       0.2222222
#> 6         0.0034246148                                       0.0000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.7111111
#> 2                                      0.2666667
#> 3                                      0.3555556
#> 4                                      0.3111111
#> 5                                      0.8444444
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                                0.2222222                   0.00000000
#> 2                                0.6666667                   0.75555556
#> 3                                0.5777778                   0.57777778
#> 4                                0.5333333                   0.88888889
#> 5                                0.8888889                   0.66666667
#> 6                                0.0000000                   0.08888889
```

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.1940665
#> 2 nonPhylo_G1_trait2                 100             0.3108555
#> 3    phylo_G1_trait1                 100             0.9748093
#> 4    phylo_G1_trait2                 100             0.9842086
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.805933481                             0
#> 2               0.670072881                             0
#> 3               0.009314317                             0
#> 4               0.004464561                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1695252
#> 2           phylo_G1_env                            0.2787971
#> 3           phylo_G1_env                            0.3838660
#> 4           phylo_G1_env                            0.3260598
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.02454136                   0.025041681
#> 2                          0.03205846                   0.035496272
#> 3                          0.59094328                   0.015734893
#> 4                          0.65814885                   0.005949253
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1       0.780891800                                            0
#> 2       0.653648193                                            0
#> 3       0.009455851                                            0
#> 4       0.009842140                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100       0.761780431
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100      -0.071022997
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100      -0.005413472
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       0.271162837
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       0.076093279
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100       0.923916821
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.05372657                  0.708053859
#> 2              -0.07466198                  0.003638987
#> 3              -0.01872596                  0.013312492
#> 4               0.28854708                 -0.017384243
#> 5               0.07291982                  0.003173455
#> 6               0.91376483                  0.010151992
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.22222222
#> 2                0.71111111                       0.57777778
#> 3                0.97777778                       0.84444444
#> 4                0.04444444                       0.04444444
#> 5                0.53333333                       0.53333333
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                            0.0000000           phylo_G1_env
#> 2                            0.6666667           phylo_G1_env
#> 3                            0.5333333           phylo_G1_env
#> 4                            0.5777778           phylo_G1_env
#> 5                            0.8888889           phylo_G1_env
#> 6                            0.0000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.04914878
#> 2                              0.01935421
#> 3                              0.06463409
#> 4                              0.12048005
#> 5                              0.10480187
#> 6                              0.27282310
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.004577788                     0.0240248096
#> 2                           -0.094016191                    -0.0005461307
#> 3                           -0.083360050                     0.0064731785
#> 4                            0.168067029                    -0.0180825270
#> 5                           -0.031882041                     0.0006665353
#> 6                            0.640941731                     0.0067273769
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1         0.6840290489                                       0.2222222
#> 2         0.0041851178                                       0.7111111
#> 3         0.0068393133                                       0.2666667
#> 4         0.0006982836                                       0.0000000
#> 5         0.0025069194                                       0.2222222
#> 6         0.0034246148                                       0.0000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                      0.7111111
#> 2                                      0.2666667
#> 3                                      0.3555556
#> 4                                      0.3111111
#> 5                                      0.8444444
#> 6                                      0.0000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                                0.2222222                   0.00000000
#> 2                                0.6666667                   0.75555556
#> 3                                0.5777778                   0.57777778
#> 4                                0.5333333                   0.88888889
#> 5                                0.8888889                   0.66666667
#> 6                                0.0000000                   0.08888889
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
#> 1    phylo_G1_env                  80             0.9871573
#> 2 phylo_G1_trait1                  80             0.9888308
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1                0.01284272                             0
#> 2                0.01057484                             0
#> [1] "Covariance results:"
#>        trait_1         trait_2 number_observations total_correlation
#> 1 phylo_G1_env phylo_G1_trait1                  80         0.7878006
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1                0.7767174                   0.01108328
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                         0                                0
#>   p_value_non_phylogenetic_correlation
#> 1                                    0
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.127524
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.126758
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.157123
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA  NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.1179
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.146984
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.121308
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.142522
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.166863
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.159641
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.154917
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.160582
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.110846
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.12771
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.128504
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.175357
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.144203
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.123194
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.093806
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.167889
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.115058
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.122907
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.139527
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.13073
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.148998
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.169702
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.093393
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.158217
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.149191
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.112859
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.152496
#>                                                                                                                                    Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, phylo_G1_env
```

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1       3.5950643
#> 2    t10       1.3589125
#> 3   t100       5.3216593
#> 4    t11       0.1944281
#> 5    t12      -3.2008335
#> 6    t13      -0.1135998
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
#> -5.5176 -0.7757  0.6811  0.9351  3.0485  5.5809      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -5.5176 -0.7021  0.6917  0.8989  2.8458  5.5809
```
