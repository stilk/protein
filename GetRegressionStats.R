GetPackages = function() {
    library(lme4)
    library(glmnet)
}


# Does mixed effect regression of Values by LogScore (i.e. mutational load) while
# controlling for cancer type and tumor purity. 
DoMixedEffectRegression = function(df) {
    model = lmer(Value ~ LogScore + (1 | type) + (1 |purity), data=df)
    ModelCoef = data.frame(summary(model)[10]$coefficients)
    ModelCoef$InterceptType = row.names(ModelCoef)
    return(ModelCoef)
}

DoRegressionPerGene = function(df, RegressionType) {
    RegressionResults = data.frame()
    AllGenes = unique(df$GeneName)
    for (i in 1:length(AllGenes)) {
        SingleGene = subset(df, df$GeneName == AllGenes[i])
        tryCatch( 
        {  if (RegressionType == 'MixedEffect') {
                RegressionResults = rbind(RegressionResults, DoMixedEffectRegression(SingleGene))
        }}, 
            error=function(error_message) { 
            message(error_message)
            print(paste0("Couldn't calculate regression statistics for ", AllGenes[i]))
        })
    }
    RegressionResults = subset(RegressionResults, RegressionResults$InterceptType != '(Intercept)')
    return(RegressionResults)
}


# DoRegressionWithinCancerType = function(df) {


# }
