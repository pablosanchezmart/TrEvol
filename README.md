
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
#> 1     t1      -3.4060594      -1.4275127   -4.8028653       0.7027674
#> 2    t10      -0.2283885       2.2471305   -0.1392432       1.8653069
#> 3   t100       2.2006552       2.6842138   -0.2853409       1.7286380
#> 4    t11       4.3962584       5.6759204    2.2349492       3.0939941
#> 5    t12      -4.6114539      -2.6385365   -2.3110741       5.5067096
#> 6    t13       0.4475639       0.4149191   -1.8555274       1.0401921
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1     -0.80670302   -1.0978337          1.1253322         1.37086171
#> 2      2.65208221    0.9924078         -0.1291769        -0.84012865
#> 3      0.09453883    2.6103445          0.6568225         1.15378084
#> 4      2.04724319    4.9572813         -0.9454392        -0.10833380
#> 5      4.12135051    3.9156076         -0.9381351        -0.48320377
#> 6      0.75295071    1.3089902         -0.2988659        -0.07249626
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1     -0.18207767         0.87567239          0.4574627       0.6722337
#> 2     -0.50458089        -1.29702694         -1.7681719      -0.4787524
#> 3      0.03392745        -0.54116560         -0.5191725      -0.4854289
#> 4     -0.06901570        -0.03181935          0.1092340      -0.1599506
#> 5     -0.74431907        -1.16954925         -0.9166246      -0.8196003
#> 6     -1.01786639         0.49377758          0.3976660       0.1262330
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

``` r
variance_covariance_results.list <- computeVarianceCovariancePartition(
  traits =  c("nonPhylo_G1_trait1", "nonPhylo_G1_trait2", "phylo_G1_trait1", "phylo_G1_trait2"),
  dataset = simulated.data$data,
  terminal_taxa = "animal",
  phylogeny = simulated.data$phylogeny
  )
```

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.3521362
#> 2 nonPhylo_G1_trait2                 100             0.2275022
#> 3    phylo_G1_trait1                 100             0.9875495
#> 4    phylo_G1_trait2                 100             0.9853158
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.647863754                             0
#> 2               0.739390260                             0
#> 3               0.006375580                             0
#> 4              -0.004504909                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100       0.871812806
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100      -0.171417837
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100      -0.002904222
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100      -0.139167146
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100      -0.078557962
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100       0.898055565
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1              0.218481987                  0.653330820
#> 2             -0.174947574                  0.003529737
#> 3              0.008426584                 -0.011330806
#> 4             -0.129640007                 -0.009527139
#> 5             -0.054938299                 -0.023619663
#> 6              0.887092986                  0.010962579
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                        0.0000000
#> 2                0.31111111                        0.3111111
#> 3                0.88888889                        0.8000000
#> 4                0.08888889                        0.1333333
#> 5                0.44444444                        0.5333333
#> 6                0.00000000                        0.0000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.7555556
#> 3                            0.6666667
#> 4                            0.8444444
#> 5                            0.3111111
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
    #> 3     ---     4   0.89

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
    #> 1     ---     2   0.65

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
```

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2745738
#> 2 nonPhylo_G1_trait2                 100             0.2943689
#> 3    phylo_G1_trait1                 100             0.9462078
#> 4    phylo_G1_trait2                 100             0.9327834
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.725426188                             0
#> 2               0.676264365                             0
#> 3               0.039031885                             0
#> 4              -0.006999637                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.2467874
#> 2           phylo_G1_env                            0.2554337
#> 3           phylo_G1_env                            0.3192331
#> 4           phylo_G1_env                            0.3132989
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.02778641                    0.03950668
#> 2                          0.03893517                    0.03277152
#> 3                          0.62697467                    0.03370815
#> 4                          0.61948449                    0.05099774
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.68591951                                            0
#> 2        0.67285963                                            0
#> 3        0.02008403                                            0
#> 4        0.01621889                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.84474110
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.07665715
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.09253999
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.24249357
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.11284430
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.89212672
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.19292469                  0.651816405
#> 2              -0.07956835                  0.002911196
#> 3              -0.07779788                 -0.014742114
#> 4              -0.23434010                 -0.008153467
#> 5              -0.10348265                 -0.009361648
#> 6               0.82225607                  0.069870651
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                        0.0000000
#> 2                0.57777778                        0.5777778
#> 3                0.62222222                        0.6222222
#> 4                0.04444444                        0.1333333
#> 5                0.13333333                        0.1333333
#> 6                0.00000000                        0.0000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                            0.0000000           phylo_G1_env
#> 2                            0.7555556           phylo_G1_env
#> 3                            0.5333333           phylo_G1_env
#> 4                            0.5333333           phylo_G1_env
#> 5                            0.5777778           phylo_G1_env
#> 6                            0.0000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                             0.162758604
#> 2                             0.035628574
#> 3                             0.003648819
#> 4                            -0.063453167
#> 5                             0.003756292
#> 6                             0.242704338
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                             0.03016609                      0.035353741
#> 2                            -0.11519692                     -0.006658836
#> 3                            -0.08144669                     -0.028951222
#> 4                            -0.17088693                     -0.003869897
#> 5                            -0.10723895                     -0.005939629
#> 6                             0.57955173                      0.057139906
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.616462664                                      0.04444444
#> 2          0.009570032                                      0.66666667
#> 3          0.014209108                                      0.97777778
#> 4         -0.004283570                                      0.17777778
#> 5         -0.003422020                                      0.97777778
#> 6          0.012730745                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                     0.26666667
#> 2                                     0.31111111
#> 3                                     0.53333333
#> 4                                     0.26666667
#> 5                                     0.04444444
#> 6                                     0.00000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                               0.08888889                    0.0000000
#> 2                               0.93333333                    0.7555556
#> 3                               0.08888889                    0.4444444
#> 4                               0.57777778                    0.7555556
#> 5                               0.44444444                    0.8888889
#> 6                               0.00000000                    0.0000000
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
    #> 3     ---     4   0.58

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
```

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1      -3.4060594
#> 2    t10      -0.2283885
#> 3   t100       2.2006552
#> 4    t11       4.3962584
#> 5    t12      -4.6114539
#> 6    t13       0.4475639
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
#> -5.9687 -1.5896  0.4348  0.1789  1.7165  5.5234      20
```

After imputation:

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -5.9687 -1.0653  0.4822  0.2785  1.7643  5.5234
```
