
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
    elif DataType == 'Protein':
        Values = GetProteinExpressionData(Dataset)
    elif DataType == 'AS':
        Values = GetAlternativeSplicingData()
    Muts = AnnotateMutationalLoad(GetMutations(Dataset), MutType=MutType)
    #Muts['NormalizedMutLoad'] = (Muts['MutLoad'] - Muts['MutLoad'].mean()) / Muts['MutLoad'].std()
    Purity = GetTumorPurity()
    Tissue = GetTissueType(Dataset)
    if DataType == 'AS':
        Muts['Barcode'] = Muts['Barcode'].str[0:12]
    Stats = Values.merge(Muts, left_on='Barcode', right_on='Barcode')
    Stats['LogScore'] = np.log10(Stats['MutLoad'] + 1)
    if Dataset == 'TCGA':
        Covariates = Tissue.merge(Purity, left_on='Barcode', right_on='Barcode')
        if DataType == 'AS':
            Covariates['Barcode'] = Covariates['Barcode'].str[0:12]
    else:
        Covariates = Tissue
    return(Stats.merge(Covariates, left_on='Barcode', right_on='Barcode'))


def GetRegressionEstimates(Dataset, DataType, MutType='SNV', GeneSet=[], ByGene=True, NormalizeY=True):
    SetUpPlottingPackages()
    df = GetRegressionStatsInput(Dataset, DataType, MutType)
    if len(GeneSet) != 0: # Subset to genes of interest
        df = df[df['GeneName'].isin(GeneSet)]
    R = ConvertPandasDFtoR(df)
    if ByGene:
        if DataType == 'Protein' or DataType == 'AS':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', False)))
        elif Dataset == 'CCLE':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', True)))
        elif Dataset == 'TCGA':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'MixedEffect', True)))
        elif Dataset == 'GTEX':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS')))
    else: # Do regression by gene group
        if Dataset == 'TCGA':
            return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGroup(R)))



def GetDeltaPSIForEachGene():
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')

Low = GetAlternativeSplicingData(AS='RI', BarcodesOfInt=Samples[Samples['MutLoad'] <= 10]['Barcode'].tolist())
High = GetAlternativeSplicingData(AS='RI', BarcodesOfInt=Samples[Samples['MutLoad'] >= 1000]['Barcode'].tolist())

# def doTTest(row):
# 	'''This is a two-sided test for the null hypothesis that 2 independent samples have identical average (expected) values.'''
# 	low = row[list(set(row.index) & set(getBarcodesOfInterest('low_SNV')))].dropna()
# 	high = row[list(set(row.index) & set(getBarcodesOfInterest('high_SNV')))].dropna()
# 	return( pd.Series({'pVal': stats.ttest_ind(low, high)[1], 
# 						'low_mean':low.mean(), 
# 						'high_mean': high.mean(), 
# 						'lowToHighDelta': low.mean() - high.mean()
# 						})
# 			)

# def calculateDeltaPSI(splice):
# 	splice = addPctMissing(splice)
# 	splice = splice[splice['pct_missing'] < 0.25]
# 	splice = splice.set_index(['symbol', 'as_id', 'splice_type', 'exons', 'from_exon', 'to_exon', 'pct_missing'])
# 	out = splice.apply(doTTest, axis=1)
# 	return(out)



# # GetNumberOfMutsAndCancerType(Dataset='TCGA')

# GetRegressionEstimates(Dataset='TCGA', DataType='AS', MutType='KsKa', GeneSet=[], ByGene=True).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ASMixedEffectRegressionEstimatesKsKaTCGA')


# GetRegressionEstimates(Dataset='TCGA', DataType='Expression', MutType='KsKa', GeneSet=[], ByGene=True).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')

# GetRegressionEstimates(Dataset='CCLE', DataType='Expression', MutType='KsKa', GeneSet=[], ByGene=True, NormalizeY=True).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ExpressionOLSRegressionEstimatesKsKaCCLE')

GetRegressionEstimates(Dataset='CCLE', DataType='Protein', MutType='KsKa', GeneSet=[], ByGene=True, NormalizeY=False).to_csv(
    '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ProteinOLSRegressionEstimatesKsKaCCLE')


# GetRegressionEstimates(Dataset='CCLE', DataType='RNAi', MutType='KsKa', GeneSet=[], ByGene=True).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/RNAiOLSRegressionEstimatesKsKaCCLE')


#RunVEPToGetDeleteriousScores(Dataset='CCLE').to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Raw/Mutations/CCLE/CCLE_Mutations_Polyphen_Annotated', sep='\t')



#GetHighlyExpressedGenes().to_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/HighlyExpressed_GTEX_GreaterThan500TPM.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='SNV').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_SNVs.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='Nonsynonymous').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_Nonsynonymous.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='Polyphen').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_Polyphen.txt')

