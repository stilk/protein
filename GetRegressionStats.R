
library(lme4)
library(glmnet)
library(lmerTest)
library(dplyr)


NormalizeValues = function(df) {
    # Z scores expression values of each gene so that they're normally distributed
    df$NormalizedValue = (df$Value - mean(df$Value))/sd(df$Value) 
	return(df)
} 


DoMixedEffectRegression = function(df, NormalizeY, DoModelDiagnostics=FALSE) { # Used for expression (TCGA)
    # Does mixed effect regression of Values by LogScore (i.e. mutational load) while
    # controlling for cancer type and tumor purity. 
    if (NormalizeY) {
         model = lmer(NormalizedValue ~ LogScore + (1 | type) + purity , data=df)
    } else {
         model = lmer(Value ~ LogScore + (1 | type) + purity, data=df)
    }
    if (DoModelDiagnostics){
        return(DoModelDiagnostics(df, model))
    } else {
        return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
    }
}


DoOLSRegressionWithGeneName = function(df, NormalizeY, DoModelDiagnostics=FALSE) { # Used for shRNA
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = data.frame(df %>% group_by(GeneName) %>% do(data.frame(NormalizeValues(.))))    
        }
        model = lm(NormalizedValue ~ LogScore + GeneName, data=df)
    } else  {
        model = lm(Value ~ LogScore + GeneName, data=df)
    }
    if (DoModelDiagnostics){
        return(DoModelDiagnostics(df, model))
    } else {
        f = summary(model)$fstatistic
        pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
        return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
    }
}

DoLinearRegression = function(df, NormalizeY, DoModelDiagnostics=FALSE) { # Used for drug and exp(CCLE)
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = NormalizeValues(df)
        }
        df = NormalizeValues(df)
        model = lm(NormalizedValue ~  LogScore, data=df)
    } else {
        model = lm(Value ~ LogScore, data=df)
    }
    if (DoModelDiagnostics){
        return(DoModelDiagnostics(df, model))
    } else {
        f = summary(model)$fstatistic
        pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
        return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
    }
}


DoRegressionPerGene = function(df, RegressionType, NormalizeY, DoModelDiagnostics=FALSE) {
    RegressionResults = data.frame()
    AllGenes = unique(df$GeneName)
    for (i in 1:length(AllGenes)) {
        SingleGene = subset(df, df$GeneName == AllGenes[i])  
        SingleGene = NormalizeValues(SingleGene)
        GeneName = unique(SingleGene$GeneName)
        tryCatch( 
        {  if (RegressionType == 'MixedEffect') {
                RegressionResults = rbind(RegressionResults, data.frame(DoMixedEffectRegression(SingleGene, NormalizeY, DoModelDiagnostics), GeneName))
        } else if (RegressionType == 'OLS') {
                RegressionResults = rbind(RegressionResults, data.frame(DoLinearRegression(SingleGene, NormalizeY, DoModelDiagnostics), GeneName))
        }}, 
            error=function(error_message) { 
            message(error_message)
            print(paste0("Couldn't calculate regression statistics for ", AllGenes[i]))
        })
    }
    return(RegressionResults)
}



DoModelDiagnostics = function(df, model) {
    # Compares fitted values vs residuals in model to assess heteroscedasticity
    Diagnostics = data.frame() # Empty df where results are appended
    df = na.omit(df) # Make sure no NAs exist, which are already omitted during regression 
    Diagnostics = rbind(Diagnostics, data.frame('PearsonsR' = cor(fitted(model),resid(model), method=c("pearson")),
        'Diagnostic' = c('FullModel'), 'Type' = c('All'))) # Compare residuals in the entire model
    Diagnostics = rbind(Diagnostics, data.frame('PearsonsR' = cor(df$LogScore, as.data.frame(summary(model)$residuals)[,1], 
        method=c("pearson")),'Diagnostic' = c('MutLoad'), 'Type' = c('All'))) # Compare residuals vs all explanatory variables
    if("purity" %in% colnames(df)) { # Only look at purity as explanatory variable in TCGA (not CCLE)
        Diagnostics = rbind(Diagnostics, data.frame('PearsonsR' = cor(df$purity, as.data.frame(summary(model)$residuals)[,1],  
            method=c("pearson")), 'Diagnostic' = c('Purity'), 'Type' = c('All')))     
        # Compare fitted vs residuals within each cancer type (i.e. full model)
        WithinGroupComp = data.frame(type = as.character(df$type), fitted = c(fitted(model)), resid = c(resid(model)) )
        WithinGroupComp = WithinGroupComp %>% group_by(type) %>% summarise(PearsonsR = cor(fitted, resid), method=c("pearson"))
        Diagnostics = rbind(Diagnostics, data.frame('PearsonsR' = WithinGroupComp$PearsonsR, 'Diagnostic' = c('FullModel'),
             'Type' = WithinGroupComp$type))   # Type also only for CCLE regression  
    }
    return(Diagnostics)
}