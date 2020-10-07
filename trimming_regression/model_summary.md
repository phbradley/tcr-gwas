# __Model Summary__: <br> Estimating effects of SNPs on V/D/J gene insertion and trimming conditional on gene

***
## Definitions and notation: 
***
* __Patient__ $i$ such that $i=1,\ldots, n$
* __T-cell__ $k$ such that $k = 1,\ldots, K_n$
* __SNP__ $j$ such that $j=1,\ldots,p$
* __Minor allele frequency__ of SNP $j$ in patient $i$ is given by $x_{ij}$
* __Gene type__ $t$ is the gene of the TCRB that we wish to condition on depending on model (either V-gene, D-gene, or J-gene).
* __Gene__ for gene_type $t$, T-cell $k$, and patient $i$ is given by $z_{t,i,k}$. 
  Each $z_{t,i,k}$ is coded as an integer such that $z_{t,i,k} = 1 \dots M_t$ where $M_t$ is the total number of candidate TCRB genes for the given gene type (either V-gene, D-gene, or J-gene)
* __Trimming or insertion length__ for T-cell $k$ and patient $i$ (on V-gene, D-gene, or J-gene depending on model) is given by $y_{i,k}$.

***
## Assumptions about the biology:
***
* Each T-cell $k$ is considered either productive (in-frame) or not-productive (out of frame).
* We measure the minor allele frequency of a SNP as a genotype with integer valued 0, 1, 2, or 3 for each individual. 
    * 0: homozygous for the dominant allele
    * 1: heterozygous for the minor allele (SNP)
    * 2: homozygous for the minor allele (SNP)
    * 3: genotype is missing for the individual
* We know the gene $z_{t,i,k}$ for gene type $t$, T-cell $k$, from patient $i$ with no error. 
* For the d-gene type, sometimes the gene call $z_{t,i,k}$ is missing (likely due to excessive trimming).
* SNPs can influence which gene $z_{t,i,k}$ is selected.
* Gene choice $z_{t,i,k}$ influences trimming and/or insertion length $y_{i,k}$.
* We are interested in the direct effect between minor allele frequency of SNP $x_{ij}$ and trimming (or insertion) length $y_{i,k}$ (removing the effect mediated by gene choice).
* An unequal number of t cells may be measured for each patient. 
* Each patient may have a patient-specific random effect on trimming (or insertion) length $y_{i,k}$ (i.e. caused by some epigenetic factor like methylation which is patient dependent, etc.)

***
## Dataset: 
***
* For the trimming analysis, for each gene type $t$ we create a "condensed dataset."
* For each patient $i \in [1,n]$, for the specified trimming type (i.e. V-gene trimming, D-gene d0 trimming, etc.), we condense the dataset by averaging the trimming by each gene $z \in [1, M_t]$ and productivity status (i.e. productive (in-frame) or not-productive (out of frame)) for the specified
  gene type $t$.
  We encode these data as a vector of length $2*M_t$ for each patient $i$. 
    * For example, if we are interested in V-gene trimming, we would condense the dataset by averaging the V-gene trimming by each possible V-gene $z \in [1, M_t]$ and productivity status (i.e. productive (in-frame) or not-productive (out of frame)).
* Since each patient likely has a different number of T-cells measured, we can let $N_{im}$ be the number of T-cells for patient $i$ that have gene $z_{t,i,k} = m$. 
  Thus, we can calculate a the proportion of T-cells with gene $z_{t,i,k} = m$ for patient $i$ as 
  $$\text{weighted gene count} = \frac{N_{im}}{\sum_{l=1}^{M_t} N_{il}}.$$
  We will use this value in the model as a weight for each observation. 
* For each patient $i$ and SNP $j$, we have a measure of minor allele frequency given by $x_{ij}$ which we can include for each entry in this condensed dataset. 
* We can concatenate these condensed datasets for all patients $i \in [1, n]$ to create the observation list (to be used in the regression) of length $2 * M_t * n$. 
* For the insertion analysis, we go through a similar data condensing procedure, except, now, instead of condensing by gene type $t$ (i.e. V, D, or J) we condense by each gene type combination (VJ, DJ, and VD) to create a "condensed dataset."
* *__Note__: for the d-gene, sometimes the gene call is missing.*
  *In these cases, for now, we exclude these entries from the observation list entirely.*
* __Example:__ Condensed dataset for V-gene and V-gene trimming. 
  We can then filter for `productive == True` or `productive == False` depending on the regression we are interested in. 
  After filtering, we will have a list of observations of length $M_t *n$ (where each row is an observation for the regression).
  
| local ID | SNP | v_gene      | productive | v_trim     | weighted_v_gene_count |
|----------|-----|-------------|------------|------------|-----------------------|
| HIP00110 | 1   | TRBV1*01    | False      | 2.9189     | 0.0017                |
| HIP00110 | 1   | TRBV1*01    | True       | 3.0112     | 0.0006                |
| HIP00110 | 1   | TRBV10-1*01 | False      | 5.4196     | 0.0025                |
| HIP00110 | 1   | TRBV10-1*01 | True       | 3.8795     | 0.0033                |
| ...      | ... | ...         | ...        | ...        | ...                   |
| HIP19717 | 2   | TRBV9*03    | False      | 1.4875     | 0.0007                |
| HIP19717 | 2   | TRBV9*03    | True       | 0.5013     | 0.0033                |

<br>
<br>

***
## Model Summary:
***
* We create a unique model for each SNP $j \in [1, p]$ and productivity status (i.e. productive (in-frame) or not-productive (out of frame)).    
* If genotypes are missing for all patients for a SNP, $j$, we do not create a model for SNP $j$.
* If all patients have the same genotype for a SNP $j\in[1,p]$, we do not create a model for SNP $j$ (since, in this case, the SNP would not have an effect on V/D/J gene insertion and trimming since the SNP genotype is the same for all observations and, thus, the slope of the regresion would be zero).
* If a patient $i \in [1,n]$ has a missing genotype for SNP $j \in [1,p],$ the patient is excluded from the observation list.
* We want to fit the model: 
  $$ y_{ik} = \beta_0 + \beta_{1j} x_{ij} + \sum_{m=1}^{M_t} \gamma_{mj} I(z_{t,i,k} = m) + \delta_i + \epsilon_{ijk} $$
  where $y_{ik}$ is the trimming length for patient $i$ and T-cell $k$, $\beta_0 + \gamma_{mj} + \delta_i$ is the intercept, (or average trimming length for T-cell with gene $z_{t,i,k} = m$, and a genotype $x_{ij}$ of 0 for SNP $j$ in patient $i$), and $\beta_{1j}$ is the slope (or direct effect of SNP $j$ on trimming length $y_{ik}$). 
* For the model, recall that since each patient likely has a different number of T-cells measured, we define a weight for each observation as $$\text{weighted gene count} = \frac{N_{im}}{\sum_{l = 1}^M N_{il}}.$$
* With this weighting and model in mind, we want to solve the following weighted least squares problem for each SNP $j$ and productivity combination: 
  $$ (\beta_0, \beta_{1j}, \gamma_{\dot j}) = \text{argmin}_{\beta_0, \beta_{1j}, \gamma_{\dot j}} \sum_{i=1}^n \sum_{m=1}^M \Bigg(\frac{N_{im}}{\sum_{l = 1}^M N_{il}}\Bigg)\big(y_{im} -(\beta_0 + \gamma_{mj} + \beta_{1j}x_{ij} + \delta_i)\big)^2 $$
* We can solve this weighted least squares problem in R using the ```lme4``` package. 
  For example, for V-gene trimming, we can model this as follows where ```localID``` is a unique identifier for each patient and ```(1|localID)``` is the term that allows for a patient specific random effect (or patient varying intercept) (Note:  an example of `condensed_data` was shown earlier in the dataset section): 
```r
regression = lmer(formula = v_trim ~ snp + v_gene + (1|localID), data = condensed_data, weights = weighted_v_gene_count)
``` 
* Using the earlier terminology, this regression will give us an $\text{intercept} = \beta_{0} + \gamma_j + \delta$ and a $\text{slope}_j = \beta_{1j}$. 
The $\text{slope}_j$ represents the direct effect of SNP $j$ on trimming length (which is what we are interested in).
* As an initial screen, using the standard error of the slope term calculated by ```lme4``` we can calculate a z-score for each SNP $j$, called $Z_j$, where 
    $$Z_j = \frac{\text{slope}_j}{SE_j}$$ 
    and compare $Z_j$ to a $N(0,1)$ distribution to get a p-value.
* If this initial screen yields a p-value less than a cutoff (which we can define), we can bootstrap across patients to get a more accurate standard error of the $\text{slope}_j$ value. 
  The code for the bootstrap is appended at the end of this write-up.
* To start, we can bootstrap 100 times to calculate the standard error, $SE_j$, of the $\text{slope}_j$, for each SNP $j$. 
    From here, as in the initial screen, we can calculate a p-value. 
* So far, we are running this protocol for a subset of SNPs which were identified in Phil's previous analyses. 
* Despite running this protocol on a select subset of SNPs, because this subset originates from an original set of 35,481,497 SNPs, we define a Bonferroni significance cutoff equal to $\frac{0.05}{35,481,497} = 1.409e-09$. 
* If after the bootstrap, the calculated p_value is less than this Bonferroni significance cutoff, we repeat the bootstrap for 1000 times to get at a more accurate $SE_j$ value, and thus, a more accurate p-value. 
* *We can find whether snps are correlated with one another by looking at their correlation coefficients based on linkage disequilibrium.*
  *If snps are adjacent with on another, and their correlation coefficients are above some cutoff which we define (for the whole block of snps in the proposed snp group), then we define this group to be a correlated cluster of snps.*


***
## Questions:
***
1. Is gene choice independent of adjacent gene choice? (i.e is V-gene choice independent of D-gene and J-gene choice?)
   * If these are not independent, should we condition on triplets in the regression (i.e. sets of VDJ genes) instead of what we are doing now (which is conditioning on only the V-gene when assessing V-gene trimming, etc.)?
2. Is it "okay" to exclude a patient from the observation list when the patient has a missing genotype for the SNP of interest (bolded bullet #4 under model summary)?
3. (Maybe an obvious biology question) How is VJ insertion possible if V and J genes aren't adjacent?

***
## Code snippets:
***

### Function for bootstrap:
Here, we assign the ```cluster_variable``` to be the patients present within the ```data``` for the indicated snp and ```trim_type``` (i.e. v gene trimming, d0 trimming, etc.).

We can repeat this protocol for the indicated number of ```repetitions```.


```r
clusboot_lmer <- function(regression, data, cluster_variable, trim_type, repetitions){
    # forms cluster variable (here, extract patient names)
    clusters <- names(table(cluster_variable))

    # Set empty matrix to hold results
    regression_results <- matrix(NA, nrow=repetitions, ncol=2)

    for(i in 1:repetitions){
        # Sample clusters (patients) with replacement
        index <- sample(1:length(clusters), length(clusters), replace=TRUE)
        # Assign clusters to index
        clusters_by_index <- clusters[index]
        # Create a contingency table for the clusters
        contingency_table_clusters <- table(clusters_by_index)
        bootdat <- NULL

        # For each cluster group (based on how many repeats in contingency table)
        for(cluster_repeat_count in 1:max(contingency_table_clusters)){
            # subset data to include only the cluster variables indicated in cluster_repeat_count
            data_subset <- data[cluster_variable %in% names(contingency_table_clusters[contingency_table_clusters %in% cluster_repeat_count]),]
            # Add data to bootdata table
            for(k in 1:cluster_repeat_count){
                bootdat <- rbind(bootdat, data_subset)
            }
        }

        # Repeat regression for this data subset, calculate standard errors
        formula = formula(regression) 

        # Set weight term for the regression
        weight = paste0("weighted_", str_split(trim_type, "_")[[1]][1], "_gene_count")

        regression_temp = lmer(formula, bootdat, weights = get(weight))
        regression_results[i,] <- c(colMeans(coef((regression_temp))$localID[1]), colMeans(coef(regression_temp)$localID[2]))
        
        # combine results into output
        coefficients_se <- cbind(c(coef(summary(regression))[1], coef(summary(regression))[2]),apply(regression_results,2,sd))
        colnames(coefficients_se) <- c("Intercept","Slope")
        rownames(coefficients_se) <- c("Estimate","Std. Error")
    } 
    return(coefficients_se)
}
```
