
import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from GetData import * 

def SetUpPlottingPackages():
    ggplot = importr('ggplot2')
    #scales = importr('scales')
    #ggpmisc = importr('ggpmisc')
    #cowplot = importr('cowplot')
    lme4 = importr('lme4')
    glmnet = importr('glmnet')
    lmertest = importr('lmerTest')
    data_table = importr('data.table')
    dplyr = importr('dplyr')
    ro.r.source('/labs/ccurtis2/tilk/scripts/protein/GetRegressionStats.R') # Source R script for plotting


def ConvertPandasDFtoR(df):
    with localconverter(ro.default_converter + pandas2ri.converter): 
        dfInR = ro.conversion.py2rpy(df) # Convert pandas dataframe to R 
    return(dfInR)


def ConvertRDataframetoPandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        dfInPd = ro.conversion.rpy2py(df) # Convert R dataframe back to Pandas
    return(dfInPd)


def GetRegressionStatsInput(Dataset, DataType, MutType):
    ''' 
    Calculates mixed effects regression (controlling for cancer type) coefficients for
    the effect of expression/drug response/rnai response to mutational load.
    @Dataset = `CCLE`, `GTEX` or `TCGA`
    @DataType = `Expression`, `RNAi` or `Drug`
    '''
    if Dataset == "TS":
        return(GetExpressionData(Dataset='TS'))
    if DataType == 'RNAi':
        Values = GetRNAiResponseData()
    elif DataType == 'Drug':
        Values = GetDrugResponseData(Screen='primary')
    elif DataType == 'Expression':
        Values = GetExpressionData(Dataset)
    Muts = AnnotateMutationalLoad(GetMutations(Dataset), MutType=MutType)
    Muts['NormalizedMutLoad'] = (Muts['MutLoad'] - Muts['MutLoad'].mean()) / Muts['MutLoad'].std()
    Purity = GetTumorPurity()
    Tissue = GetTissueType(Dataset)
    Stats = Values.merge(Muts, left_on='Barcode', right_on='Barcode')
    #Stats['LogScore'] = np.log10(Stats['MutLoad'] + 1)
    if Dataset == 'TCGA':
        Covariates = Tissue.merge(Purity, left_on='Barcode', right_on='Barcode')
    else:
        Covariates = Tissue
    return(Stats.merge(Covariates, left_on='Barcode', right_on='Barcode'))


def GetRegressionEstimates(Dataset, DataType, MutType='SNV', GeneSet=[], ByGene=True):
    SetUpPlottingPackages()
    df = GetRegressionStatsInput(Dataset, DataType, MutType)
    if len(GeneSet) != 0: # Subset to genes of interest
        df = df[df['GeneName'].isin(GeneSet)]
    R = ConvertPandasDFtoR(df)
    if ByGene:
        if Dataset == 'TCGA':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'MixedEffect')))
        elif Dataset == 'CCLE':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS')))
        elif Dataset == 'GTEX':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS')))
    else: # Do regression by gene group
        if Dataset == 'TCGA':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGroup(R)))

def GetRegressionToGeneSetsOfInterest(Group, Dataset, DataType, MutType):
    Dir = os.getcwd() + '/GeneSets/'
    RegressionResults = pd.DataFrame()
    if Group == 'Chaperome':
        set = pd.read_csv(Dir + 'Human_Chaperome_TableS1A_PMID_29293508', sep='\t').rename(columns={'Gene':'GeneName','Level2':'Group'})
    for GeneGroup in set['Group'].unique():
        GeneSet = set[set['Group'] == GeneGroup]['GeneName']
        RegressionResults = RegressionResults.append(GetRegressionEstimates(Dataset, DataType, MutType, GeneSet, ByGene=False).assign(GeneGroup = GeneGroup))
    return(RegressionResults)



GetRegressionEstimates(Dataset='TCGA', DataType='Expression', MutType='SNV', GeneSet=[], ByGene=True).to_csv(
    '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/StandardizedExpressionMixedEffectRegressionEstimatesTCGA')



#GetHighlyExpressedGenes().to_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/HighlyExpressed_GTEX_GreaterThan500TPM.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='SNV').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_SNVs.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='Nonsynonymous').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_Nonsynonymous.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='Polyphen').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_Polyphen.txt')

