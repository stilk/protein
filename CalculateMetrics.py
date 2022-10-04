'''
This script uses the rpy2 package to perform regression analyses (on all different data types and datasets) in python using code in R. 
'''


import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import os.path
from GetRCodeIntoPython import *
from GetAnnotations import *
from GetData import * 




def GetRegressionStatsInput(Dataset, DataType, MutType, AllDrugs=False):
    ''' 
    Outputs a dataframe used for regression analysis which compares how a variable of interest (e.g. expression values) 
    changes in response to to mutational load. Functions for aggregating input data are imported from @GetData.
    @Dataset = `CCLE` or `TCGA`
    @DataType = `Expression`, `RNAi`, `Drug`, `AS` (Alternative Splicing) or `Protein`
    @AllDrugs = Optional variable for drug input data, if false will only report proteostasis drug
    '''
    if DataType == 'RNAi':
        Values = GetRNAiResponseData()
    elif DataType == 'Drug':
        Values = GetDrugResponseData('primary', AllDrugs)
    elif DataType == 'Expression':
        Values = GetExpressionData(Dataset)
    elif DataType == 'Protein':
        Values = GetProteinExpressionData(Dataset)
    elif DataType == 'AS':
        Values = GetAlternativeSplicingData()
    Muts = AnnotateMutationalLoad(GetPointMutations(Dataset), MutType=MutType)
    Purity = GetTumorPurity(Dataset)
    Tissue = GetTissueType(Dataset)
    if DataType == 'AS':
        Muts['Barcode'] = Muts['Barcode'].str[0:12]
    Stats = Values.merge(Muts, left_on='Barcode', right_on='Barcode')
    Stats['LogScore'] = np.log10(Stats['MutLoad'] + 1)
    if Dataset == 'TCGA':
        Covariates = Tissue.merge(Purity, left_on='Barcode', right_on='Barcode')
        Covariates['ShortBarcode'] = Covariates['Barcode'].str[0:12]
        Covariates = Covariates.merge(GetPatientAge(), left_on='ShortBarcode', right_on='ShortBarcode', how='left')
    else:
        Covariates = Tissue
    return(Stats.merge(Covariates, left_on='Barcode', right_on='Barcode'))



def GetProteinRegression(Dataset, DataType='Protein', MutType='KsKa', ModelDiagnostics='None', ShuffleTMB=False):  
    '''
    Uses linear regression to measure the association of expression and mutational load. Input data for regressinon is 
    taken from @GetRegressionInput. 
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @ModelDiagnostics = If 'SummarizedByAllGenes', will compare fitted values vs residuals in a given regression model to assess heteroscedasticity.
                        summarized as correlation coef. If 'RawValsByComplex', will output raw residuals for genes in complexes of interest.
                        Model diagnostics are performed by @GetRegressionModelDiagnostics. If 'None', will not perform model diagnostics.
    '''
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset, DataType, MutType)
    if ShuffleTMB: # Used to compare distribution of null vs empirical p-values for empirical FDR in Fig.1?
        df['LogScore'] = df.groupby(['type','GeneName'])['LogScore'].transform(np.random.permutation) # shuffle
    R = ConvertPandasDFtoR(df)
    if Dataset == 'CCLE':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', True, ModelDiagnostics)))
    elif Dataset == 'TCGA':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLSWithPurity', True, ModelDiagnostics)))
    elif Dataset == 'CPTAC':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'MixedWithoutPurity', True, ModelDiagnostics)))


def GetExpressionRegression(Dataset, DataType='Expression', MutType='KsKa', ModelDiagnostics='None', ShuffleTMB=False):  
    '''
    Uses linear regression to measure the association of expression and mutational load. Input data for regressinon is 
    taken from @GetRegressionInput. Depending on the dataset, will use OLS or GLMM (mixed effect) regression.
    Outputs dataframe of regression results.
    @Dataset = `CCLE` or `TCGA` 
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @ModelDiagnostics = If 'SummarizedByAllGenes', will compare fitted values vs residuals in a given regression model to assess heteroscedasticity.
                        summarized as correlation coef. If 'RawValsByComplex', will output raw residuals for genes in complexes of interest.
                        Model diagnostics are performed by @GetRegressionModelDiagnostics. If 'None', will not perform model diagnostics.
    '''
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset, DataType, MutType)
    if ShuffleTMB: # Used to compare distribution of null vs empirical p-values for empirical FDR in Fig.1?
        df['LogScore'] = df.groupby(['type','GeneName'])['LogScore'].transform(np.random.permutation) # shuffle
    R = ConvertPandasDFtoR(df)
    if Dataset == 'CCLE':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', True, ModelDiagnostics)))
    elif Dataset == 'TCGA':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'MixedEffect', True, ModelDiagnostics)))
    

def DoRegressionByAge(DataType='Expression', MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/'):
    '''
    Uses function @GetExpressionRegression to perform linear regression to measure the association of 
    expression and mutational load by each age group in TCGA. Outputs dataframe of regression results.
    @DataType = `Expression`,`RNAi`, `Drug`, `AS` (Alternative Splicing) or `Protein`
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @OutDir = Output directory to write to file
    '''
    Output = pd.DataFrame() # Where results are appended to
    SetUpRegressionPackages()
    df = GetRegressionStatsInput('TCGA', DataType, MutType)
    df = df.merge( GetGeneAnnotationsOfInterest(), left_on='GeneName', right_on='Hugo')
    df['AgeBin'] = pd.cut(df['age'], [9,30,40,50,60,70,80,90], labels=['10-30','30-40','40-50','50-60','60-70','70-80','80-90'])
    for AgeGroup in df['AgeBin'].unique():
        AgeSubset = df[df['AgeBin'] == AgeGroup]
        RegResult = ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(ConvertPandasDFtoR(AgeSubset), 'MixedEffect', True, False))
        RegResult['AgeBin'] = AgeGroup
        Output = Output.append(RegResult)
    return(Output)
    Output.to_csv(OutDir + 'TCGA_GLMM_ByAgeGroups')


def GetshRNARegression(ModelDiagnostics='None'):
    '''
    Uses linear regression to measure the association of viability after expression knockdown and mutational load. 
    Input data for regressinon is taken from @GetRegressionInput. Outputs dataframe of regression results.
    @Dataset = `CCLE` or `TCGA` 
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @ModelDiagnostics = If 'SummarizedByAllGenes', will compare fitted values vs residuals in a given regression model to assess heteroscedasticity.
                        summarized as correlation coef. If 'RawValsByComplex', will output raw residuals for genes in complexes of interest.
                        Model diagnostics are performed by @GetRegressionModelDiagnostics. If 'None', will not perform model diagnostics.
    '''
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset='CCLE', DataType='RNAi', MutType='KsKa')
    groups = GetGeneAnnotationsOfInterest()
    df = df.merge(groups, left_on='GeneName', right_on='Hugo')
    Out = pd.DataFrame()
    for Group in groups['Group'].unique():
        SingleGroup = ConvertPandasDFtoR(df[df['Group'] == Group].dropna()[['LogScore','Value','GeneName']])
        Reg = ConvertRDataframetoPandas(ro.r.DoOLSRegressionWithGeneName(SingleGroup, False, ModelDiagnostics)).reset_index()
        Reg['Group'] = Group
        Reg['subgroup'] = groups[groups['Group'] == Group]['subgroup'].unique()[0]
        Out = Out.append(Reg)
    Out = Out.dropna().rename(columns={'index':'Coefficient'})
    if ModelDiagnostics != 'None':
        Out = Out[Out['Coefficient'] == 'LogScore']
    return(Out)


def DoPurityRegression(ModelDiagnostics='None'): 
    ''' Reviewer request to do GLMM by purity smdfh'''
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset='TCGA', DataType='Expression', MutType='KsKa', AllDrugs=False)
    R = ConvertPandasDFtoR(df)
    return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'ByPurity', True, 'None')).reset_index())
    

def GetPerDrugRegression(ModelDiagnostics='None'):
    '''
    Uses linear regression to measure the association of viability after drug inhibition and mutational load. 
    Input data for regressinon is taken from @GetRegressionInput. Outputs dataframe of regression results.
    @Dataset = `CCLE` or `TCGA` 
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @ModelDiagnostics = If 'SummarizedByAllGenes', will compare fitted values vs residuals in a given regression model to assess heteroscedasticity.
                        summarized as correlation coef. If 'RawValsByComplex', will output raw residuals for genes in complexes of interest.
                        Model diagnostics are performed by @GetRegressionModelDiagnostics. If 'None', will not perform model diagnostics.
    '''
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset='CCLE', DataType='Drug', MutType='KsKa')
    Out = pd.DataFrame()
    for Drug in df['name'].unique():
        SingleDrug = ConvertPandasDFtoR(df[df['name'] == Drug].dropna())
        Reg= ConvertRDataframetoPandas(ro.r.DoLinearRegression(SingleDrug, True, ModelDiagnostics)).reset_index()
        Reg = pd.concat([df.query('name == @Drug')[['name','subgroup','Group','screenType']].head(2).reset_index(), 
            Reg.rename(columns={'index':'Coefficient'})], axis=1) ## Add drug name and group to regression output
        Out = Out.append(Reg)
    return(Out)


def GetWithinCancerType(Dataset, DataType, MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/WithinCancerType/'):
    ''' Need to add functionality for different DataTypes? '''
    SetUpRegressionPackages()
    AllCancerTypes = GetRegressionStatsInput(Dataset, DataType, MutType) # Get input data for 
    for CancerType in AllCancerTypes['type'].unique(): 
        Knifed = AllCancerTypes[AllCancerTypes['type'] == CancerType] # Select one cancer type from loop
        R = ConvertPandasDFtoR(Knifed) # Convert df to R for regression
        print('Analyzing cancer type... ' + str(CancerType))
        if Dataset == 'CCLE':
            OutFile = OutDir + Dataset + DataType + MutType + 'OLSRegressionWithinCancerType' + CancerType.replace(' ','_').replace('/','_') 
            if DataType == 'RNAi':
                GetshRNARegression(Knifed).assign(CancerType = CancerType).to_csv(OutFile)
            elif DataType == 'Drug': # Do grouped regression 
                GetPerDrugRegression(Knifed).assign(CancerType = CancerType).to_csv(OutFile)
            else: # Expression or protein, which we're performing per gene regression on
                OneTypeKnifed = ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', True)).assign(
                    CancerType = CancerType).to_csv(OutFile)
        elif Dataset == 'TCGA':
            OutFile = OutDir + Dataset + DataType + MutType + 'MixedEffectRegressionWithinCancerType' + CancerType
            OneTypeKnifed = ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLSWithPurity', True)).assign(
                CancerType = CancerType).to_csv(OutFile)


def GetWithinCancerTypeGrouped(Dataset, DataType, MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/WithinCancerType/'):
    ''' Adds a grouped gene name for OLS  '''
    SetUpRegressionPackages()
    groups = GetGeneAnnotationsOfInterest() # Get groups grouped regression will be done on
    AllCancerTypes = GetRegressionStatsInput(Dataset, DataType, MutType) # Get input data for regression
    for CancerType in AllCancerTypes['type'].unique().tolist() + ['All Cancers']: # Loop through individual cancer types and all cancer types as ref point
        OutFile = OutDir + Dataset + DataType + MutType + 'ByComplexGroupRegressionWithinCancerTypeIncludingGeneNameAndPurity' + CancerType.replace(' ','_').replace('/','_') 
        if CancerType != 'All Cancers':
            Knifed = AllCancerTypes[AllCancerTypes['type'] == CancerType] 
        else: # Don't need to subset if looking at all cancer types
            Knifed = AllCancerTypes
        Knifed = Knifed.merge(groups, left_on='GeneName', right_on='Hugo')
        R = ConvertPandasDFtoR(Knifed) # Convert df to R for regression
        print('Analyzing cancer type... ' + str(CancerType)) 
        Out = pd.DataFrame()
        for Group in groups['Group'].unique():
            SingleGroup = ConvertPandasDFtoR(Knifed[Knifed['Group'] == Group].dropna()[['LogScore','Value','GeneName','purity','type']])
            if Dataset == 'TCGA':
                if CancerType != 'All Cancers':
                    Reg = ConvertRDataframetoPandas(ro.r.DoLinearRegressionWithPurityAndGeneName(SingleGroup, True, ModelDiagnostics='None')).reset_index()
                else:
                    Reg = ConvertRDataframetoPandas(ro.r.DoMixedEffectRegression(SingleGroup, True, ModelDiagnostics='None')).reset_index()
            elif Dataset == 'CCLE':
                Reg = ConvertRDataframetoPandas(ro.r.DoOLSRegressionWithGeneName(SingleGroup, True, ModelDiagnostics='None')).reset_index()
            Reg['Group'] = Group
            Reg['subgroup'] = groups[groups['Group'] == Group]['subgroup'].unique()[0]
            Out = Out.append(Reg)
        Out = Out.dropna().rename(columns={'index':'Coefficient'})
        Out.to_csv(OutFile)
        print('Done with cancer type: ' + CancerType)
    

def JacknifeAcrossCancerTypes(Dataset, DataType, CancerTypeToRemove='', MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/Jacknife/'):
    ''' 
    Removes a cancer type from the dataset and re-calculates regression coefficients.
    Writes jacknifed regression coefficients to file since files per gene are huge. Only set up for per gene for now. 
    @Dataset = `CCLE` or `TCGA` 
    @DataType = `Expression`,`RNAi`, `Drug`, `AS` (Alternative Splicing) or `Protein`
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @CancerTypeToRemove = a string of which individual cancer type to remove. If empty loops through
        all cancer types. If not, removes cancer types of that string and exists loop after 1 iteration.
    @OutDir = Output directory to write to file
    '''
    SetUpRegressionPackages()
    AllCancerTypes = GetRegressionStatsInput(Dataset, DataType, MutType) # Get input data for all cancer types
    for CancerType in AllCancerTypes['type'].unique(): # Loop through all cancer types
        if CancerTypeToRemove != '': # Don't loop through all cancer types; only do jacknife for one type
            CancerType = CancerTypeToRemove
        Knifed = AllCancerTypes[AllCancerTypes['type'] != CancerType] # Remove cancer type from loop
        R = ConvertPandasDFtoR(Knifed) # Convert df to R for regression
        print('Removing cancer type... ' + str(CancerType))
        if Dataset == 'CCLE':
            OutFile = OutDir + Dataset + DataType + MutType + 'OLSRegressionJacknifed' + CancerType.replace(' ','_').replace('/','_') 
            if os.path.isfile(OutFile): # If file already exists, skip
                continue
            if DataType == 'RNAi':
                GetshRNARegression(Knifed).assign(CancerTypeRemoved = CancerType).to_csv(OutFile)
            elif DataType == 'Drug': # Do grouped regression 
                GetPerDrugRegression(Knifed).assign(CancerTypeRemoved = CancerType).to_csv(OutFile)
            else: # Expression or protein, which we're performing per gene regression on
                OneTypeKnifed = ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', True)).assign(
                    CancerTypeRemoved = CancerType).to_csv(OutFile)
        elif Dataset == 'TCGA':
            OutFile = OutDir + Dataset + DataType + MutType + 'MixedEffectRegressionJacknifed' + CancerType
            if os.path.isfile(OutFile): # If file already exists, skip
                continue
            OneTypeKnifed = ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'MixedEffect', True)).assign(
                CancerTypeRemoved = CancerType).to_csv(OutFile)
        if CancerTypeToRemove != '': # Only iterating on one cancer type; break loop if true
            print('Done with analysis for cancer type: ' + str(CancerType))
            break


def GetRegressionModelDiagnostics(Dataset, DataType, MutType='KsKa', ModelDiagnostics='SummarizedByAllGenes', OutDir=os.getcwd() + '/Data/Regression/Diagnostics/'):
    '''
    Calculates Pearson's R for fitted values with residuals for the full regression model, 
    and residuals vs all explanatory variables. Outputs statistics for each gene. Outputs dataframe of diagnostics results.
    @Dataset = `CCLE` or `TCGA` 
    @DataType = `Expression`,`RNAi`, `Drug`, `AS` (Alternative Splicing) or `Protein`
    @MutTtype = variable denoting how to define mutational load (default is all protein coding mutations, i.e. 'KsKa')
    @AllGenes = string when 'SummarizedByAllGenes' will generate correlation coeff of residuals vs fitted values for all genes, 
                if 'RawValsByComplex' will output raw residuals for genes in complexes of interest.
    @OutDir = Output directory to write to file
    '''
    print(Dataset, DataType, MutType,ModelDiagnostics)
    if DataType == 'Expression':
        df = GetExpressionRegression(Dataset, DataType='Expression', MutType=MutType, ModelDiagnostics = ModelDiagnostics)
    elif DataType == 'shRNA':
        df = GetshRNARegression(ModelDiagnostics = ModelDiagnostics)
    elif DataType == 'Drug':
        df = GetPerDrugRegression(ModelDiagnostics = ModelDiagnostics)
    df.to_csv(OutDir + Dataset + DataType + MutType + 'ModelDiagnosticResidualVsFittedVals' + ModelDiagnostics)
    return(df)
 




### PROTEIN ###

# GetProteinRegression().to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/MedianProteinOLSRegressionEstimatesKsKaCCLE')
# GetProteinRegression('CPTAC').to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ProteinMixedEffectRegressionWithoutPurityEstimatesKsKaCPTAC')




# ran within cancer type expression
# run jacknife protein

# GetWithinCancerType(Dataset='TCGA', DataType='Expression')


DoPurityRegression().to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/TCGA_ByPurity_MixedEffect_RR')


# GetWithinCancerType(Dataset='CCLE', DataType='Protein')
# GetWithinCancerTypeGrouped('TCGA','Expression')
  
# JacknifeAcrossCancerTypes('CCLE', 'Protein')

# import random
# RandStr = str(random.randint(0, 900000000))
# print("STARTING NOW")
# print(RandStr)

# GetExpressionRegression(Dataset='TCGA', DataType='Expression', MutType='KsKa', ModelDiagnostics='None', ShuffleTMB=True).to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/Shuffling/ExpressionMixedEffectRegressionEstimatesTCGAShuffledTMB_' + RandStr)


# GetRegressionModelDiagnostics(Dataset='TCGA', DataType='Expression', MutType='KsKa', ModelDiagnostics='Homoscedasticity').to_csv('/labs/ccurtis2/tilk/scripts/protein/homoscesdascity-TCGA')


# GetRegressionModelDiagnostics(Dataset='TCGA', DataType='Expression', MutType='KsKa', ModelDiagnostics='DHARMA').to_csv('/labs/ccurtis2/tilk/scripts/protein/DHARMA-TCGA')

#GetRegressionModelDiagnostics(Dataset='TCGA', DataType='Expression', MutType='KsKa', ModelDiagnostics='WithinCancerType').to_csv('/labs/ccurtis2/tilk/scripts/protein/Skew-WithinCancerType-TCGA')



# GetRegressionModelDiagnostics(Dataset='TCGA', DataType='Expression', MutType='KsKa', ModelDiagnostics='NoAutocorrelation').to_csv('/labs/ccurtis2/tilk/scripts/protein/autocor-TCGA')


