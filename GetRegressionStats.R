GetPackages = function() {
    library(lme4)
    library(glmnet)
    library(lmerTest)
}

NormalizeValues = function(df) {
    # Z scores expression values of each gene so that they're normally distributed
    df$NormalizedValue = (df$Value - mean(df$Value))/sd(df$Value) 
	return(df)
} 


DoMixedEffectRegression = function(df) {
    # Does mixed effect regression of Values by LogScore (i.e. mutational load) while
    # controlling for cancer type and tumor purity. 
    model = lmer(NormalizedValue ~ LogScore + (1 | type) + (1 |purity), data=df)
    return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
}

DoLinearRegression = function(df) {
    model = lmer(Value ~ LogScore, data=df)
    f = summary(model)$fstatistic
    pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
    return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
}


DoRegressionPerGene = function(df, RegressionType) {
    RegressionResults = data.frame()
    AllGenes = unique(df$GeneName)
    for (i in 1:length(AllGenes)) {
        SingleGene = subset(df, df$GeneName == AllGenes[i])  
        SingleGene = NormalizeValues(SingleGene)
        GeneName = unique(SingleGene$GeneName)
        tryCatch( 
        {  if (RegressionType == 'MixedEffect') {
                RegressionResults = rbind(RegressionResults, data.frame(DoMixedEffectRegression(SingleGene), GeneName))
        } else if (RegressionType == 'OLS') {
                RegressionResults = rbind(RegressionResults, data.frame(DoLinearRegression(SingleGene), GeneName))
        }}, 
            error=function(error_message) { 
            message(error_message)
            print(paste0("Couldn't calculate regression statistics for ", AllGenes[i]))
        })
    }
    return(RegressionResults)
}


# DoRegressionWithinCancerType = function(df) {


# }
