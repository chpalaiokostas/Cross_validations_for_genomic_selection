# Wrap up function for cross validations in genomic selection models using BGLR

-----------------------------------------
# The following arguments can be supplied
# A data frame with phenotypes: Y
# The first column should be the trait of interest
# The rest of the columns should be the fixed effect
# The function expects each of fixed effect to be already of the appropriate type
# E.g for factors change the variables appropriately beforehand 
# Number of cross validations (default set to 5): cv
# Number of iterations for the McMC (default set to 1e06): numItter 
# The first 10% of the iteratins is considered as burn in  
# Genotypes in the form of matrix. No missing values are allowed.
# Instead of genotypes a relationship matrix can be provided. 
# Type of model: "BRR", "BayesA", "BayesB", "BayesC", "BL","RKHS"
# More info about the above models can be found in BGLR documentation
# The type of trait: trait_type
# The allowed values for trait type are either "gaussian" or "ordinal"
# Whether the trait is censored: censored
# The allowed values for censored are TRUE or FALSE
# Whether a relationship matrix is provided instead of a matrix with genotypes 
# The allowed values for relationship_matrix are TRUE or FALSE   

-----------------------------------------
# Returns the Pearson correlation between 
# estimated breeding values and the phenotypes from each CV set
# If the trait is categorical probabilities are returned  
# Also returns the bias of the estimated breeding values from each CV set 
# This is estimated by regressing the adjusted for fixed effects phenotype
# on the estimated breeding value  
  
cross_validation <- function(Y,
                             cv=5,numItter=1e06,
                             modelMatrix,
                             model,trait_type,
                             censored=F,relationship_matrix=F) 
  {
    N <- nrow(Y)
    y <- Y[,1]
    fixed_effects <- gsub("\n","",as.character(formula(Y))[3])
    adj_pheno <- resid(lm(as.formula(paste(y,"~",fixed_effects))))
    group_size <- floor(N/fold)
    ncv <- cumsum(rep(group_size,cv))
    
    result <- c()
    if(relationship_matrix) {
      ETA <- list(FIXED=list(as.formula(paste("~",fixed_effects)),
                             model="FIXED", data=Y),
                  MRK=list(K=modelMatrix,model=model))  
    }
    else {
      ETA <- list(FIXED=list(as.formula(paste("~",fixed_effects)),
                             model="FIXED", data=Y),
                  MRK=list(X=modelMatrix,model=model)) 
    }
    #start_index <- 0
    for(i in 1:cv) {
      rem <- 1:N
      for(j in 1:length(ncv)) {
        yNA <- y
        tst <- sample(rem,size=group_size,replace=F)
        rem <- setdiff(rem,tst)
        yNA[tst] <- NA
        if(censored) {
          yNA <- log(yNA)
          a <- rep(NA, length(y))
          b <- rep(NA, length(y))
          a[index] <- 3.9
          b[index] <- Inf
          yNA[index] <- NA
          a[tst] <- -Inf
          b[tst] <- Inf
          fm<- BGLR(y=yNA,response_type=trait_type,a=a,b=b,
                    ETA=ETA,nIter=numItter,burnIn=0.1*numItter,thin=100)
          a <- rep(NA, length(y))
          b <- rep(NA, length(y))
          reg_coef <- coef(lm(adj_pheno[tst]) ~ fm$yHat[tst])
        }
        else {
          fm <- BGLR(y=yNA,response_type=trait_type,ETA=ETA,
                     nIter=numItter,burnIn=0.1*numItter,thin=100)
          reg_coef <- coef(lm(adj_pheno[tst] ~ fm$yHat[tst]))
        }
        if(trait_type=="ordinal") {
          temp <- cbind(fm$probs[tst,2],y[tst])
          result <- rbind(result,temp)
        }
        else {
          reg_coef <- coef(lm(adj_pheno[tst] ~ fm$yHat[tst]))
          result <- c(result,cor(fm$yHat[tst],y[tst]),reg_coef)
        }
      }
    }
    result
}

