
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
#> 1     t1        3.998936      4.23434847     3.278168      -0.3786966
#> 2    t10        2.487570      0.08915715     1.347419       6.1368644
#> 3   t100        6.918799      4.13158618     5.064814       1.1295081
#> 4    t11       -1.837807     -3.08484866    -2.940020       1.5416422
#> 5    t12       -5.728792     -5.09234043    -3.043322       4.9911319
#> 6    t13        5.372559      2.81630876     5.130768       1.4225013
#>   phylo_G2_trait2 phylo_G2_env nonPhylo_G1_trait1 nonPhylo_G1_trait2
#> 1      -1.4597684   -2.0448966         -0.9713297         -0.7264323
#> 2       5.1130869    5.6779172         -0.2493384         -0.3796915
#> 3      -0.5122594    2.3182360         -1.5362784         -1.2641065
#> 4       1.9090098   -1.7396322          0.9720901          1.7560260
#> 5       4.0524621    2.2465310         -0.5040611         -0.7329585
#> 6      -0.2075244    0.8932458          0.7694495          2.1038458
#>   nonPhylo_G1_env nonPhylo_G2_trait1 nonPhylo_G2_trait2 nonPhylo_G2_env
#> 1      -1.5978500          0.7884529         0.95635936       0.5882107
#> 2       0.4959047         -0.1339331         0.14344631       1.0531433
#> 3      -0.7700887         -0.2759235        -0.63602537       0.6724618
#> 4       1.6581907          0.6041272         0.64245481       1.2908262
#> 5      -0.2678424         -0.1140744         0.01380234      -0.7724040
#> 6       1.6970502          0.4091042         1.41277877       0.9997337
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
    #> 1 nonPhylo_G1_trait1                 100             0.2005389
    #> 2 nonPhylo_G1_trait2                 100             0.2153287
    #> 3    phylo_G1_trait1                 100             0.9952457
    #> 4    phylo_G1_trait2                 100             0.9933262
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.799461126                             0
    #> 2               0.759853900                             0
    #> 3               0.007661004                             0
    #> 4               0.009820479                             0
    #> [1] "Covariance results:"
    #>              trait_1            trait_2 number_observations total_correlation
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100      0.8468717178
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100     -0.1132058477
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100     -0.0007881875
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100     -0.0517228826
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100      0.0335425394
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100      0.8961409325
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.07788367                  0.768988046
    #> 2              -0.14238503                  0.029179180
    #> 3              -0.02865554                  0.027867357
    #> 4              -0.06626352                  0.014540637
    #> 5               0.01580965                  0.017732886
    #> 6               0.89163557                  0.004505361
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                 0.0000000                       0.04444444
    #> 2                 0.2666667                       0.13333333
    #> 3                 0.9777778                       0.93333333
    #> 4                 0.6666667                       0.57777778
    #> 5                 0.7111111                       0.71111111
    #> 6                 0.0000000                       0.00000000
    #>   p_value_non_phylogenetic_correlation
    #> 1                            0.0000000
    #> 2                            0.0000000
    #> 3                            0.1333333
    #> 4                            0.3555556
    #> 5                            0.2222222
    #> 6                            0.0000000

We can now explore the variance-covariance results. First, let’s look at
the variance results.

``` r
variance_covariance_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.2005389
#> 2 nonPhylo_G1_trait2                 100             0.2153287
#> 3    phylo_G1_trait1                 100             0.9952457
#> 4    phylo_G1_trait2                 100             0.9933262
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1               0.799461126                             0
#> 2               0.759853900                             0
#> 3               0.007661004                             0
#> 4               0.009820479                             0
```

Now, let’s look at the covariance results, which are reported as
correlation.

``` r
variance_covariance_results.list$covarianceResults
#>              trait_1            trait_2 number_observations total_correlation
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100      0.8468717178
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100     -0.1132058477
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100     -0.0007881875
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100     -0.0517228826
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100      0.0335425394
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100      0.8961409325
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.07788367                  0.768988046
#> 2              -0.14238503                  0.029179180
#> 3              -0.02865554                  0.027867357
#> 4              -0.06626352                  0.014540637
#> 5               0.01580965                  0.017732886
#> 6               0.89163557                  0.004505361
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                 0.0000000                       0.04444444
#> 2                 0.2666667                       0.13333333
#> 3                 0.9777778                       0.93333333
#> 4                 0.6666667                       0.57777778
#> 5                 0.7111111                       0.71111111
#> 6                 0.0000000                       0.00000000
#>   p_value_non_phylogenetic_correlation
#> 1                            0.0000000
#> 2                            0.0000000
#> 3                            0.1333333
#> 4                            0.3555556
#> 5                            0.2222222
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
            correlation_threshold = 0.3)
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
    #> 1     ---     2   0.77

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
    #> 1 nonPhylo_G1_trait1                 100             0.1613379
    #> 2 nonPhylo_G1_trait2                 100             0.2144486
    #> 3    phylo_G1_trait1                 100             0.9885536
    #> 4    phylo_G1_trait2                 100             0.9823994
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1                0.83866211                             0
    #> 2                0.76203352                             0
    #> 3                0.00647835                             0
    #> 4                0.01532113                             0
    #>   environmental_variable non_attributed_phylogenetic_variance
    #> 1           phylo_G1_env                            0.1551215
    #> 2           phylo_G1_env                            0.1651586
    #> 3           phylo_G1_env                            0.3108448
    #> 4           phylo_G1_env                            0.3833945
    #>   environmental_phylogenetic_variance labile_environmental_variance
    #> 1                         0.006216392                   0.041542393
    #> 2                         0.049289966                   0.068415554
    #> 3                         0.677708777                   0.002812152
    #> 4                         0.599004976                   0.006219307
    #>   residual_variance p_value_non_attributed_phylogenetic_variance
    #> 1       0.797119714                                            0
    #> 2       0.717135862                                            0
    #> 3       0.008634228                                            0
    #> 4       0.011381251                                            0
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
    #> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.82704619
    #> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.10704399
    #> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.06952799
    #> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.01646330
    #> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.10386288
    #> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90582927
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1               0.05887379                  0.768172407
    #> 2               0.10477135                  0.002272641
    #> 3               0.03278277                  0.036745221
    #> 4              -0.02323361                  0.006770309
    #> 5              -0.12150601                  0.017643124
    #> 6               0.89937532                  0.006453944
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                0.00000000                       0.04444444
    #> 2                0.26666667                       0.22222222
    #> 3                0.44444444                       0.66666667
    #> 4                0.71111111                       0.71111111
    #> 5                0.08888889                       0.04444444
    #> 6                0.00000000                       0.00000000
    #>   p_value_non_phylogenetic_correlation environmental_variable
    #> 1                           0.00000000           phylo_G1_env
    #> 2                           0.88888889           phylo_G1_env
    #> 3                           0.04444444           phylo_G1_env
    #> 4                           0.66666667           phylo_G1_env
    #> 5                           0.31111111           phylo_G1_env
    #> 6                           0.00000000           phylo_G1_env
    #>   non_attributed_phylogenetic_correlation
    #> 1                              0.05761414
    #> 2                              0.10241141
    #> 3                             -0.02656435
    #> 4                              0.08472478
    #> 5                              0.09891656
    #> 6                              0.24664880
    #>   environmental_phylogenetic_correlation labile_environmental_correlation
    #> 1                            0.001259652                      0.037818138
    #> 2                            0.002359936                      0.005915452
    #> 3                            0.059347117                      0.008068386
    #> 4                           -0.107958391                      0.008571655
    #> 5                           -0.220422567                      0.012445676
    #> 6                            0.652726522                      0.002529999
    #>   residual_correlation p_value_non_attributed_phylogenetic_correlation
    #> 1          0.730354270                                      0.08888889
    #> 2         -0.003642811                                      0.04444444
    #> 3          0.028676835                                      0.62222222
    #> 4         -0.001801346                                      0.04444444
    #> 5          0.005197448                                      0.00000000
    #> 6          0.003923945                                      0.00000000
    #>   p_value_environmental_phylogenetic_correlation
    #> 1                                     0.97777778
    #> 2                                     0.93333333
    #> 3                                     0.40000000
    #> 4                                     0.08888889
    #> 5                                     0.04444444
    #> 6                                     0.00000000
    #>   p_value_labile_environmental_correlation p_value_residual_correlation
    #> 1                               0.08888889                   0.00000000
    #> 2                               0.26666667                   0.75555556
    #> 3                               0.57777778                   0.17777778
    #> 4                               0.35555556                   0.93333333
    #> 5                               0.17777778                   0.62222222
    #> 6                               0.22222222                   0.04444444

We can see now these components in both variance and covariance results.

``` r
variance_covariance_environment_results.list$varianceResults
#>                trait number_observations phylogenetic_variance
#> 1 nonPhylo_G1_trait1                 100             0.1613379
#> 2 nonPhylo_G1_trait2                 100             0.2144486
#> 3    phylo_G1_trait1                 100             0.9885536
#> 4    phylo_G1_trait2                 100             0.9823994
#>   non_phylogenetic_variance p_value_phylogenetic_variance
#> 1                0.83866211                             0
#> 2                0.76203352                             0
#> 3                0.00647835                             0
#> 4                0.01532113                             0
#>   environmental_variable non_attributed_phylogenetic_variance
#> 1           phylo_G1_env                            0.1551215
#> 2           phylo_G1_env                            0.1651586
#> 3           phylo_G1_env                            0.3108448
#> 4           phylo_G1_env                            0.3833945
#>   environmental_phylogenetic_variance labile_environmental_variance
#> 1                         0.006216392                   0.041542393
#> 2                         0.049289966                   0.068415554
#> 3                         0.677708777                   0.002812152
#> 4                         0.599004976                   0.006219307
#>   residual_variance p_value_non_attributed_phylogenetic_variance
#> 1       0.797119714                                            0
#> 2       0.717135862                                            0
#> 3       0.008634228                                            0
#> 4       0.011381251                                            0
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
#> 1 nonPhylo_G1_trait1 nonPhylo_G1_trait2                 100        0.82704619
#> 2 nonPhylo_G1_trait1    phylo_G1_trait1                 100        0.10704399
#> 3 nonPhylo_G1_trait1    phylo_G1_trait2                 100        0.06952799
#> 4 nonPhylo_G1_trait2    phylo_G1_trait1                 100       -0.01646330
#> 5 nonPhylo_G1_trait2    phylo_G1_trait2                 100       -0.10386288
#> 6    phylo_G1_trait1    phylo_G1_trait2                 100        0.90582927
#>   phylogenetic_correlation non_phylogenetic_correlation
#> 1               0.05887379                  0.768172407
#> 2               0.10477135                  0.002272641
#> 3               0.03278277                  0.036745221
#> 4              -0.02323361                  0.006770309
#> 5              -0.12150601                  0.017643124
#> 6               0.89937532                  0.006453944
#>   p_value_total_correlation p_value_phylogenetic_correlation
#> 1                0.00000000                       0.04444444
#> 2                0.26666667                       0.22222222
#> 3                0.44444444                       0.66666667
#> 4                0.71111111                       0.71111111
#> 5                0.08888889                       0.04444444
#> 6                0.00000000                       0.00000000
#>   p_value_non_phylogenetic_correlation environmental_variable
#> 1                           0.00000000           phylo_G1_env
#> 2                           0.88888889           phylo_G1_env
#> 3                           0.04444444           phylo_G1_env
#> 4                           0.66666667           phylo_G1_env
#> 5                           0.31111111           phylo_G1_env
#> 6                           0.00000000           phylo_G1_env
#>   non_attributed_phylogenetic_correlation
#> 1                              0.05761414
#> 2                              0.10241141
#> 3                             -0.02656435
#> 4                              0.08472478
#> 5                              0.09891656
#> 6                              0.24664880
#>   environmental_phylogenetic_correlation labile_environmental_correlation
#> 1                            0.001259652                      0.037818138
#> 2                            0.002359936                      0.005915452
#> 3                            0.059347117                      0.008068386
#> 4                           -0.107958391                      0.008571655
#> 5                           -0.220422567                      0.012445676
#> 6                            0.652726522                      0.002529999
#>   residual_correlation p_value_non_attributed_phylogenetic_correlation
#> 1          0.730354270                                      0.08888889
#> 2         -0.003642811                                      0.04444444
#> 3          0.028676835                                      0.62222222
#> 4         -0.001801346                                      0.04444444
#> 5          0.005197448                                      0.00000000
#> 6          0.003923945                                      0.00000000
#>   p_value_environmental_phylogenetic_correlation
#> 1                                     0.97777778
#> 2                                     0.93333333
#> 3                                     0.40000000
#> 4                                     0.08888889
#> 5                                     0.04444444
#> 6                                     0.00000000
#>   p_value_labile_environmental_correlation p_value_residual_correlation
#> 1                               0.08888889                   0.00000000
#> 2                               0.26666667                   0.75555556
#> 3                               0.57777778                   0.17777778
#> 4                               0.35555556                   0.93333333
#> 5                               0.17777778                   0.62222222
#> 6                               0.22222222                   0.04444444
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
    #> 3     ---     4   0.65

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
    #> 1    phylo_G1_env                  80             0.9890101
    #> 2 phylo_G1_trait1                  80             0.9966911
    #>   non_phylogenetic_variance p_value_phylogenetic_variance
    #> 1               0.010989934                             0
    #> 2               0.003321337                             0
    #> [1] "Covariance results:"
    #>        trait_1         trait_2 number_observations total_correlation
    #> 1 phylo_G1_env phylo_G1_trait1                  80         0.7876999
    #>   phylogenetic_correlation non_phylogenetic_correlation
    #> 1                0.7838179                  0.003882023
    #>   p_value_total_correlation p_value_phylogenetic_correlation
    #> 1                         0                                0
    #>   p_value_non_phylogenetic_correlation
    #> 1                           0.04444444
    #> Loading required package: missForest
    #> Loading required package: doRNG
    #> Loading required package: foreach
    #> Loading required package: rngtools
    #> Loading required package: randomForest
    #> randomForest 4.7-1.1
    #> Type rfNews() to see new features/changes/bug fixes.
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.571291
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.582397
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.564929
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.585528
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.623845
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.563091
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.643653
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.564944
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.625289
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.550123
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.620964
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.570902
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.561874
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.579123
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA   NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.58888
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.583286
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.567972
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.597185
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.637971
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.625456
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.548591
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.656681
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.646684
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.631591
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.661131
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.599512
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.609132
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.584385
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.534201
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env
    #>          Variable   N N_Obs N_NA    NRMSE
    #> 1 phylo_G1_trait1 100    80   20 1.627407
    #>                                                                                                                                                                                                                                           Model
    #> 1 ( phylo_G1_trait1 ) <-  Phylo_axis_1, Phylo_axis_2, Phylo_axis_3, Phylo_axis_4, Phylo_axis_5, Phylo_axis_6, Phylo_axis_7, Phylo_axis_8, Phylo_axis_9, Phylo_axis_10, Phylo_axis_11, Phylo_axis_12, Phylo_axis_13, Phylo_axis_14, phylo_G1_env

``` r
head(imputed.data$round3$ximp)
#>   animal phylo_G1_trait1
#> 1     t1        3.998936
#> 2    t10        2.487570
#> 3   t100        6.918799
#> 4    t11       -1.837807
#> 5    t12       -5.728792
#> 6    t13        5.372559
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
#>  -8.767  -4.809  -1.413  -1.514   1.203   7.037      20
```

After imputation

``` r
summary(imputed.data$phylo_G1_trait1_imputed)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  -8.767  -4.601  -1.401  -1.399   1.008   7.037
```
