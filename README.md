Wrap up function for cross validations in genomic selection models using BGLR

-----------------------------------------
The following arguments can be supplied
* A data frame with phenotypes: Y
* The first column should be the trait of interest
* The rest of the columns should be the fixed effect
* The function expects each of fixed effect to be already of the appropriate type.
E.g for factors change the variables appropriately beforehand 
* Number of cross validations (default set to 5): cv
* Number of iterations for the McMC (default set to 1e06): numItter. 
The first 10% of the iteratins is considered as burn in  
* Genotypes in the form of matrix. No missing values are allowed.
Instead of genotypes a relationship matrix can be provided. 
* Type of model: "BRR", "BayesA", "BayesB", "BayesC", "BL","RKHS".
More info about the above models can be found in BGLR documentation
* The type of trait: trait_type.
The allowed values for trait type are either "gaussian" or "ordinal"
* Whether the trait is censored: censored.
The allowed values for censored are TRUE or FALSE
* Whether a relationship matrix is provided instead of a matrix with genotypes. 
The allowed values for relationship_matrix are TRUE or FALSE   

-----------------------------------------
Returns the following:

* Pearson correlation between estimated breeding values and the phenotypes from each CV set
* If the trait is categorical probabilities are returned  
* Also returns the bias of the estimated breeding values from each CV set. This is estimated by regressing
the adjusted for fixed effects phenotype on the estimated breeding value  
