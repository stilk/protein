
library(lme4)
library(glmnet)
library(lmerTest)
library(dplyr)


NormalizeValues = function(df) {
    # Z scores expression values of each gene so that they're normally distributed
    df$NormalizedValue = (df$Value - mean(df$Value))/sd(df$Value) 
	return(df)
} 


DoMixedEffectRegression = function(df, NormalizeY) {
    # Does mixed effect regression of Values by LogScore (i.e. mutational load) while
    # controlling for cancer type and tumor purity. 
    if (NormalizeY) {
         model = lmer(NormalizedValue ~ LogScore + (1 | type) + purity , data=df)
    } else {
         model = lmer(Value ~ LogScore + (1 | type) + purity, data=df)
    }
    return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
}

DoMixedEffectRegressionWithGeneName = function(df) {
    model = lmer(NormalizedValue ~ LogScore + (1 | type) + purity + GeneName, data=df)
    return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
}


DoOLSRegressionWithGeneName = function(df, NormalizeY) {
    if (NormalizeY) {
        model = lm(NormalizedValue ~ LogScore + GeneName, data=df)
    } else  {
        model = lm(Value ~ LogScore + GeneName, data=df)
    }
    f = summary(model)$fstatistic
    pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
    return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
}

DoLinearRegression = function(df, NormalizeY) {
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = NormalizeValues(df)
        }
        df = NormalizeValues(df)
        model = lm(NormalizedValue ~  LogScore, data=df)
    } else {
        model = lm(Value ~ LogScore, data=df)
    }
    f = summary(model)$fstatistic
    pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
    return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
}


DoRegressionPerGene = function(df, RegressionType, NormalizeY) {
    RegressionResults = data.frame()
    AllGenes = unique(df$GeneName)
    for (i in 1:length(AllGenes)) {
        SingleGene = subset(df, df$GeneName == AllGenes[i])  
        SingleGene = NormalizeValues(SingleGene)
        GeneName = unique(SingleGene$GeneName)
        tryCatch( 
        {  if (RegressionType == 'MixedEffect') {
                RegressionResults = rbind(RegressionResults, data.frame(DoMixedEffectRegression(SingleGene, NormalizeY), GeneName))
        } else if (RegressionType == 'OLS') {
                RegressionResults = rbind(RegressionResults, data.frame(DoLinearRegression(SingleGene, NormalizeY), GeneName))
        }}, 
            error=function(error_message) { 
            message(error_message)
            print(paste0("Couldn't calculate regression statistics for ", AllGenes[i]))
        })
    }
    return(RegressionResults)
}


DoRegressionPerGroup = function(df, RegressionType, NormalizeY) {
    if (NormalizeY) {
        print('wrong')
        df = data.frame(df %>% group_by(GeneName) %>% do(data.frame(NormalizeValues(.))))    
    }
    tryCatch( 
    {   if (RegressionType == 'MixedEffect') {
            RegressionResults = DoMixedEffectRegressionWithGeneName(df)
        } else if (RegressionType == 'OLS') {
            print(df)
            RegressionResults = DoOLSRegressionWithGeneName(df, NormalizeY)
        }
    }, 
        error=function(error_message) { 
        message(error_message) # Print error statement
        print(paste0("Couldn't calculate regression statistics for group."))
    })
    return(RegressionResults)
}
