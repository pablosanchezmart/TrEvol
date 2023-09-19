
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
#> 1     t1     -0.12551213       0.2744017   -1.1484200       1.1713674
#> 2    t10      1.65936416       1.7133124   -0.5936037       1.3638125
#> 3   t100      1.97169830       1.8249651   -0.6046426      -0.3355446
#> 4    t11      0.85555915       1.4615536   -3.2296305       1.5086043
#> 5    t12      4.54201374       4.6127613    3.0544885       2.4223415
#> 6    t13     -0.07377813      -1.2903494    0.4711629       0.9644056
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1        2.474330    0.3885456         -0.8608908         -1.0329551
#> 2        2.007978    2.6573042          0.4327287          0.6683132
#> 3       -1.817849   -1.9131215         -0.5306932         -1.4525249
#> 4        1.300341    0.6962146          1.1433915          0.4820141
#> 5        3.054268    3.5338250         -0.3090415         -0.1922071
#> 6        1.492739    1.8312989         -1.2726714         -0.7556224
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1       0.1806628         -0.6876218         -0.9099886      -0.2829271
#> 2      -0.3144536          1.1219661          1.0231629       1.5371690
#> 3      -1.1732056         -1.7396443         -1.5721464      -0.8184107
#> 4       1.4521710          0.1094516          1.2943483       0.8237995
#> 5       0.5010888          0.2068468          0.6668777       0.4280024
#> 6      -1.1396247          0.7656596          0.7491775       0.6025300
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
#> 1 nonPhylo_G1_trait1                 100             0.1302592
#> 2 nonPhylo_G1_trait2                 100             0.1482165
#> 3    phylo_G1_trait1                 100             0.9612777
#> 4    phylo_G1_trait2                 100             0.9705843
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.869740781                             0
#> 2               0.827367132                             0
#> 3               0.028518757                             0
#> 4               0.009874016                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.87474052
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.09070165
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.05169835
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.12473385
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100        0.03371805
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90371633
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1              0.070932511                   0.80380801
#> 2              0.068801888                   0.02189976
#> 3             -0.074235455                   0.02253711
#> 4             -0.200936209                   0.07620235
#> 5             -0.007239775                   0.04095783
#> 6              0.860988915                   0.04272742
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                        0.0000000
#> 2                0.13333333                        0.4888889
#> 3                0.66666667                        0.5777778
#> 4                0.04444444                        0.0000000
#> 5                0.57777778                        0.8888889
#> 6                0.00000000                        0.0000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.2666667
#> 3                            0.5777778
#> 4                            0.0000000
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
    #> 3     ---     4   0.86

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
    #> 1     ---     2   0.8

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
#> 1 nonPhylo_G1_trait1                 100             0.1883986
#> 2 nonPhylo_G1_trait2                 100             0.2322678
#> 3    phylo_G1_trait1                 100             0.9519148
#> 4    phylo_G1_trait2                 100             0.9619232
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1                0.81160145                             0
#> 2                0.74445388                             0
#> 3                0.04417617                             0
#> 4                0.02054654                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1364071
#> 2           phylo_G1_env                            0.2117111
#> 3           phylo_G1_env                            0.2291957
#> 4           phylo_G1_env                            0.2993947
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                          0.05199141                    0.10665431
#> 2                          0.02055672                    0.06462310
#> 3                          0.72271910                    0.02575143
#> 4                          0.66252845                    0.01980486
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1        0.70494713                                            0
#> 2        0.70310907                                            0
#> 3        0.02233378                                            0
#> 4        0.01827194                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.83287195
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100       -0.09935269
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100       -0.10820506
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.00224301
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.02733812
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.92155533
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.11580300                   0.71706894
#> 2              -0.15993517                   0.06058248
#> 3              -0.15773114                   0.04952608
#> 4              -0.02717197                   0.02492896
#> 5              -0.05768465                   0.03034653
#> 6               0.87754108                   0.04401425
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                       0.08888889
#> 2                 0.1333333                       0.04444444
#> 3                 0.1333333                       0.04444444
#> 4                 0.8888889                       0.71111111
#> 5                 0.8444444                       0.48888889
#> 6                 0.0000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                           0.00000000           phylo_G1_env
#> 2                           0.08888889           phylo_G1_env
#> 3                           0.17777778           phylo_G1_env
#> 4                           0.40000000           phylo_G1_env
#> 5                           0.44444444           phylo_G1_env
#> 6                           0.17777778           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.07092588
#> 2                             -0.02887145
#> 3                             -0.03807022
#> 4                             -0.08206280
#> 5                             -0.04125355
#> 6                              0.19989263
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                             0.04487713                      0.109113155
#> 2                            -0.13106372                      0.017936990
#> 3                            -0.11966092                      0.023012799
#> 4                             0.05489083                      0.003061639
#> 5                            -0.01643111                      0.024727806
#> 6                             0.67764845                      0.038244180
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.607955788                                      0.08888889
#> 2          0.042645489                                      0.80000000
#> 3          0.026513278                                      0.57777778
#> 4          0.021867323                                      0.00000000
#> 5          0.005618728                                      0.26666667
#> 6          0.005770066                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                     0.35555556
#> 2                                     0.08888889
#> 3                                     0.31111111
#> 4                                     0.66666667
#> 5                                     0.80000000
#> 6                                     0.00000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                               0.04444444                   0.00000000
#> 2                               0.31111111                   0.04444444
#> 3                               0.44444444                   0.08888889
#> 4                               0.88888889                   0.44444444
#> 5                               0.44444444                   0.71111111
#> 6                               0.26666667                   0.26666667
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
    #> 3     ---     4   0.68

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
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.017276
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.016571
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.016291
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA     NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.9871928
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.020098
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.039054
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.029353
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.050496
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.022244
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.042841
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA     NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.9813492
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.042293
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA     NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.9856614
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.024516
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA     NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.9943353
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.005228
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.992658
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA     NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.9985526
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.032098
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.046421
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA     NRMSE
#> 1 phylo_G1_trait1 100    80   20 0.9954985
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.016861
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.01535
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.060895
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.011157
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.036735
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.011023
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.018348
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA    NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.024584
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
#>          Variable   N N_Obs N_NA   NRMSE
#> 1 phylo_G1_trait1 100    80   20 1.02537
#>                                                                                                                                                                                                             Model
#> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, phylo_G1_env
```

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1     -0.12551213
#> 2    t10      1.65936416
#> 3   t100      1.97169830
#> 4    t11      0.85555915
#> 5    t12      4.04933634
#> 6    t13     -0.07377813
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
#> -4.7760 -1.7943 -0.4834 -0.2556  1.0099  8.0186      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -4.7760 -1.6380 -0.4753 -0.1473  1.0240  8.0186
```
