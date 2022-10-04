# This script runs all of the GLMs and GLMMs used in the manuscript.
# This output of this R code is imported into Python for downstream analysis in CalculateMetrics.py

library(lme4)
library(glmnet)
library(lmerTest)
library(dplyr)
library(moments) 
library(car)
library(nortest)
library(lmtest)

NormalizeValues = function(df) {
    # Z scores expression values in the independent variable so that they're normally distributed.
    # Adds column 'NormalizedValue' to dataframe from input column 'Value'
    # Parameters:
    #   @df = dataframe of inputs for normalization. 
    df$NormalizedValue = (df$Value - mean(df$Value))/sd(df$Value) 
	return(df)
} 


DoMixedEffectRegression = function(df, NormalizeY, ModelDiagnostics='None') { # Used for expression (TCGA)
    # Does mixed effect regression of Values by LogScore (i.e. mutational load) while controlling for cancer type and
    # tumor purity. Outputs a dataframe of regression results. 
    # Parameters:
    #   @df = dataframe of raw inputs used for regression. Should have columns for the dependent variables named
    #       'purity','type', and 'LogScore'. Independent variable column is named 'Value'.
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = data.frame(df %>% group_by(GeneName) %>% do(data.frame(NormalizeValues(.))))    
        }
        model = lmer(NormalizedValue ~ LogScore + (1 | type) + purity , data=df)
    } else {
        model = lmer(Value ~ LogScore + (1 | type) + purity, data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
    }
}


DoRegressionByPurity = function(df, NormalizeY, ModelDiagnostics='None') {
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = data.frame(df %>% group_by(GeneName) %>% do(data.frame(NormalizeValues(.))))    
        }
        model = lmer(purity ~ (1 | type) + NormalizedValue , data=df)
    } else {
        model = lmer(purity ~  (1 | type) + Value, data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
    }
}
    


DoMixedEffectRegressionWithoutPurity = function(df, NormalizeY, ModelDiagnostics='None') { # Used for expression (TCGA)
    # Does mixed effect regression of Values by LogScore (i.e. mutational load) while controlling for cancer type and
    # tumor purity. Outputs a dataframe of regression results. 
    # Parameters:
    #   @df = dataframe of raw inputs used for regression. Should have columns for the dependent variables named
    #       'purity','type', and 'LogScore'. Independent variable column is named 'Value'.
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = data.frame(df %>% group_by(GeneName) %>% do(data.frame(NormalizeValues(.))))    
        }
        model = lmer(NormalizedValue ~ LogScore + (1 | type) , data=df)
    } else {
        model = lmer(Value ~ LogScore + (1 | type) , data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        return(data.frame(summary(model)[10]$coefficients, summary(model)[11]))
    }
}


DoOLSRegressionWithGeneName = function(df, NormalizeY, ModelDiagnostics='None') { # Used for shRNA
    # Does OLS regression of 'Value' by 'LogScore' and 'GeneName'. Outputs a dataframe of regression results.
    # Parameters:
    #   @df = dataframe of raw inputs used for regression. Should have columns for the dependent variables named
    #       'GeneName', and 'LogScore'. Independent variable column is named 'Value'.
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = data.frame(df %>% group_by(GeneName) %>% do(data.frame(NormalizeValues(.))))    
        }
        model = lm(NormalizedValue ~ LogScore + GeneName, data=df)
    } else  {
        model = lm(Value ~ LogScore + GeneName, data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        f = summary(model)$fstatistic
        pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
        return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
    }
}





DoLinearRegressionWithPurityAndGeneName = function(df, NormalizeY, ModelDiagnostics='None') { # Used for within cancer type grouped TCGA
    # Does OLS regression of Values by LogScore (i.e. mutational load) and outputs a dataframe of regression results. 
    # Parameters:
    #   @df = dataframe of raw inputs used for regression. Should have column for the dependent variable named 'LogScore'. 
    #       Independent variable column is named 'Value'.
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = NormalizeValues(df)
        }
        df = NormalizeValues(df)
        model = lm(NormalizedValue ~  LogScore + purity + GeneName, data=df)
    } else { # do not normalize Y values
        model = lm(Value ~ LogScore + purity + GeneName, data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        f = summary(model)$fstatistic
        pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
        return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
    }
}


DoLinearRegressionWithPurity = function(df, NormalizeY, ModelDiagnostics='None') { # Used for protein TCGA
    # Does OLS regression of Values by LogScore (i.e. mutational load) and outputs a dataframe of regression results. 
    # Parameters:
    #   @df = dataframe of raw inputs used for regression. Should have column for the dependent variable named 'LogScore'. 
    #       Independent variable column is named 'Value'.
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = NormalizeValues(df)
        }
        df = NormalizeValues(df)
        model = lm(NormalizedValue ~  LogScore + purity, data=df)
    } else { # do not normalize Y values
        model = lm(Value ~ LogScore + purity, data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        f = summary(model)$fstatistic
        pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
        return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
    }
}


DoLinearRegression = function(df, NormalizeY, ModelDiagnostics='None') { # Used for drug and exp(CCLE)
    # Does OLS regression of Values by LogScore (i.e. mutational load) and outputs a dataframe of regression results. 
    # Parameters:
    #   @df = dataframe of raw inputs used for regression. Should have column for the dependent variable named 'LogScore'. 
    #       Independent variable column is named 'Value'.
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    if (NormalizeY) {
        if(!("NormalizedValue" %in% names(df))) { # If normalized values don't exist in dataframe, normalize vals
           df = NormalizeValues(df)
        }
        df = NormalizeValues(df)
        model = lm(NormalizedValue ~  LogScore, data=df)
    } else { # do not normalize Y values
        model = lm(Value ~ LogScore, data=df)
    }
    if (ModelDiagnostics != 'None'){ # Do model diagnostics
        return(DoModelDiagnostics(df, model,ModelDiagnostics))
    } else { # Don't do model diagnostics
        f = summary(model)$fstatistic
        pVal = pf(f[1],f[2],f[3],lower.tail=F) # p value for the overall fit of the model againt null that R2=0
        return(data.frame(summary(model)[9], summary(model)[6], summary(model)[4]$coefficients, pVal))
    }
}


DoRegressionPerGene = function(df, RegressionType, NormalizeY, ModelDiagnostics='None') {
    # Takes a dataframe of values (e.g., expression/shRNA/drug) and mutational load for each tumor/cell line and 
    # runs @DoLinearRegression of @DoMixedEffectsRegression by each gene in the dataframe.
    # Parameters:
    #   @df = dataframe of raw inputs used for regression
    #   @RegressionType = a string of which regression model to use 'MixedEffect' or 'OLS'
    #   @NormalizeY = boolean of whether to normalize y values using @NormalizeValues
    #   @ModelDiagnostics = string of which type of model diagnostic to run by @DoModelDiagnostic,
    #       will not run if string is 'None' 
    RegressionResults = data.frame()
    AllGenes = unique(df$GeneName)
    for (i in 1:length(AllGenes)) {
        SingleGene = subset(df, df$GeneName == AllGenes[i])  
        SingleGene = NormalizeValues(SingleGene)
        GeneName = unique(SingleGene$GeneName)
        tryCatch( 
        {  if (RegressionType == 'MixedEffect') {
            RegressionResults = rbind(RegressionResults, data.frame(DoMixedEffectRegression(SingleGene, NormalizeY, ModelDiagnostics), GeneName))
        } else if (RegressionType == 'OLS') {
            RegressionResults = rbind(RegressionResults, data.frame(DoLinearRegression(SingleGene, NormalizeY, ModelDiagnostics), GeneName))
        } else if (RegressionType == 'OLSWithPurity') {
            RegressionResults = rbind(RegressionResults, data.frame(DoLinearRegressionWithPurity(SingleGene, NormalizeY, ModelDiagnostics), GeneName))
        } else if (RegressionType == 'MixedWithoutPurity') {
            RegressionResults = rbind(RegressionResults, data.frame(DoMixedEffectRegressionWithoutPurity(SingleGene, NormalizeY, ModelDiagnostics), GeneName))
        } else if (RegressionType == 'ByPurity') {
            RegressionResults = rbind(RegressionResults, data.frame(DoRegressionByPurity(SingleGene, NormalizeY, ModelDiagnostics), GeneName))
        }}, 
            error=function(error_message) { 
            print(message(error_message))
            print(paste0("Couldn't calculate regression statistics for ", AllGenes[i]))
        })
    }
    return(RegressionResults)
}



DoModelDiagnostics = function(df, model, DiagnosticType='') {
    # Compares fitted values vs residuals in a given regression model to assess heteroscedasticity.
    # Compares residuals in the entire model overall, and residuls against all explanatory variables.
    # Parameters:
    #   @df = dataframe of regression model inputs used (e.g., expression values/load/purity/type for each tumor)
    #   @model = raw regression model output that contains residuals
    #   @DiagnosticType = If 'SummarizedByAllGenes', will compare fitted values vs residuals in a given regression model to assess heteroscedasticity.
    #                    summarized as correlation coef. If 'RawValsByComplex', will output raw residuals for genes in complexes of interest.
    #                   Model diagnostics are performed by @GetRegressionModelDiagnostics. If 'None', will not perform model diagnostics.

    Diagnostics = data.frame() # Empty df where results are appended
    df = na.omit(df) # Make sure no NAs exist, which are already omitted during regression 
    if (DiagnosticType == 'SummarizedByAllGenes') { # Pearson's R 
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
    } else if (DiagnosticType == 'RawValsByComplex'){
        Diagnostics = data.frame(gene = as.character(df$GeneName), type = as.character(df$type), 
            fitted = c(fitted(model)), resid = c(resid(model)) )  # Raw residuals for every gene and cancer type
    } else if (DiagnosticType == 'NoAutocorrelation') { # independence between fitted vs predicted values
        # Perform Kolmogorov-Smirnov test to see if two samples come from same distribution, would expect to significantly reject null if not correlated
        # Perform Durbin Watson test for detecting autocorrelation from residuals in regression model
        # Durbin Watson test should be a value of 2 if no autocorrelatoin
        # Should perhaps also run Breusch-Godfrey test?
        Diagnostics = data.frame(durbin_watson = durbinWatsonTest(c(resid(model))), # Durbin Watson test should be a value of 2 if no autocorrelation
            ks_statistic = ks.test(c(resid(model)), c)$statistic, ks_pvalue = ks.test(c(resid(model)), c(fitted(model)))$p.value) 
    } else if (DiagnosticType == 'Homoscedasticity') { # residuals should follow normal distribution with zero mean and equal variance
        # Perform Anderson-Darling test where null hypothesis is that it is normally distributed; expect to not reject null
        # Add measurements of skewness (test of symmetry in distr) and kurtosis (measure whether distribution is shifted; should be 3 if normally dist)
        # Can't add Breusch-Pagan Test for heteroscedasticity; only applies to OLS

        Diagnostics = data.frame(ad_statistic = ad.test((c(resid(model))))$p.value, 
                ad_pvalue=ad.test((c(resid(model))))$p.value,
                skewness = skewness(c(resid(model))), 
                kurtosis = kurtosis(c(resid(model))))
                
    } else if (DiagnosticType == 'AllModelDiagnostics') {
        df = data.frame( ks.test(resid(model), model@frame$NormalizedValue))
            #homogeneity of variance within groups 
            # VarPerGroup = data.frame(AbsSqrdError = abs(resid(model))^2, Type = model@frame$type)
            # data.frame(anova(lm(AbsSqrdError ~ Type, VarPerGroup)))[1:1,]
            # Test whether variance of the residuals is equal in groups; 
            # want to REJECT null model if assumption of homoscedasticity is met
    } else if (DiagnosticType == 'DHARMA') {
        library(DHARMa)
        simulationOutput = simulateResiduals(fittedModel = model, plot = F)
        Diagnostics = data.frame(
            #BartlettsTest=bartlett.test(formula, dataset), 
            Uniformity_PVal=testUniformity(simulationOutput)$p.value, #  tests if the overall distribution conforms to expectations; NS
            Uniformity_Statistic = testUniformity(simulationOutput)$statistic, #  tests if the overall distribution conforms to expectations; NS
            Outlier_Pval = testOutliers(simulationOutput)$p.value, # tests if there are more simulation outliers than expected
            Dispersion_PVal = testDispersion(simulationOutput)$p.value, # tests if the simulated dispersion is equal to the observed dispersion
            Dispersion_Statistic = testDispersion(simulationOutput)$statistic, # tests if the simulated dispersion is equal to the observed dispersion
            HomogeneityGroupsFValue = testCategorical(simulationOutput, catPred = model@frame$type)$homogeneity[2]$"F value"[1], #tests residuals against a categorical predictor
            HomogeneityGroupsPVal = testCategorical(simulationOutput, catPred = model@frame$type)$homogeneity[3]$"Pr(>F)"[1], #tests residuals against a categorical predictor
            testZeroInflation(simulationOutput)$p.value # tests if there are more zeros in the data than expected from the simulations
        )
    } else if (DiagnosticType == 'WithinCancerType') {
        Diagnostics = as.data.frame(data.frame(type = model@frame$type, residual = resid(model)) %>% group_by(type) %>% 
            summarize(skew= skewness(residual)))
    }
    return(Diagnostics)
}



DoGeneSetEnrichment = function(GeneSet) {
    # Outputs a dataframe of gene set enrichment analysis results
    # Parameters:
    #   @GeneSet = a vector of gene names to query
    library('gprofiler2')
    result = data.frame(gost(GeneSet)$result)
    return(result[c('source','term_name','p_value')])
}

