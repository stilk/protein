
import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from GetData import * 

def SetUpRegressionPackages():
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
        Covariates['ShortBarcode'] = Covariates['Barcode'].str[0:12]
        Covariates = Covariates.merge(GetPatientAge(), left_on='ShortBarcode', right_on='ShortBarcode', how='left')
        if DataType == 'AS':
            Covariates['Barcode'] = Covariates['Barcode'].str[0:12]
    else:
        Covariates = Tissue
    return(Stats.merge(Covariates, left_on='Barcode', right_on='Barcode'))


def GetRegressionEstimates(Dataset, DataType, MutType='SNV', GeneSet=[], ByGene=True, NormalizeY=True):
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset, DataType, MutType)
    if len(GeneSet) != 0: # Subset to genes of interest
        df = df[df['GeneName'].isin(GeneSet)]
    R = ConvertPandasDFtoR(df)
    if ByGene:
        if DataType == 'Protein' or DataType == 'AS' or DataType == 'RNAi':
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
        #elif Dataset == 'CCLE':


def GetSubsetOfGeneAnnotationsOfInterest():
    corum = corum = pd.read_csv('/labs/ccurtis2/tilk/02_DNDS/separateDatabaseFiles/CORUM/coreComplexes.txt.zip', sep='\t')
    chap = pd.read_csv('/labs/ccurtis2/tilk/09_PROTEIN/Human_Chaperome_TableS1A_PMID_29293508', sep='\t')
    groups = pd.concat([
        pd.DataFrame({'Group': 'Cytoplasmic Ribosomes','subgroup' : 'Translation', 'Hugo' : corum[corum['ComplexName'].str.contains('subunit, cytoplasmic')]['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].drop_duplicates().tolist()}),
        pd.DataFrame({'Group': 'Mitochondrial Ribosomes', 'subgroup' : 'Translation', 'Hugo' : corum[corum['ComplexName'].str.contains('subunit, mitochondrial')]['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].drop_duplicates().tolist()}),
        pd.DataFrame({'Group': 'Mitochondrial Chaperones' , 'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'MITO') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': '20S Core' , 'subgroup' : 'Proteasome', 'Hugo' : corum[corum['ComplexName'] == '20S proteasome']['subunits(Gene name)'].str.split(';', expand=True).melt()['value']  }),
        pd.DataFrame({'Group': '19S Regulatory Particle' , 'subgroup' : 'Proteasome', 'Hugo' : ['PSMC1','PSMC2','PSMC3','PSMC4','PSMC5','PSMC6', 'PSMD1','PSMD2','PSMD3','PSMD4','PSMD6','PSMD7','PSMD8','PSMD11','PSMD12','PSMD13','PSMD14'] }),
        pd.DataFrame({'Group': 'Small HS' , 'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'sHSP') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 100' ,'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'HSP100') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 90', 'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'HSP90') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 70','subgroup' : 'Chaperones', 'Hugo' : chap[(chap['Level2'] == 'HSP70') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 60','subgroup' : 'Chaperones', 'Hugo' : chap[(chap['Level2'] == 'HSP60') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 40', 'subgroup' : 'Chaperones','Hugo' : chap[(chap['Level2'] == 'HSP40') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']})])
    return(groups)

def GetshRNARegression():
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset='CCLE', DataType='RNAi', MutType='KsKa')
    groups = GetSubsetOfGeneAnnotationsOfInterest()
    df = df.merge(groups, left_on='GeneName', right_on='Hugo')
    Out = pd.DataFrame()
    for Group in groups['Group'].unique():
        SingleGroup = ConvertPandasDFtoR(df[df['Group'] == Group].dropna()[['LogScore','Value','GeneName']])
        #Reg= ConvertRDataframetoPandas(ro.r.DoLinearRegression(SingleGroup, NormalizeY=True)).reset_index()
        Reg = ConvertRDataframetoPandas(ro.r.DoOLSRegressionWithGeneName(SingleGroup, False)).reset_index()
        Reg = pd.concat([df.query('Group == @Group')[['subgroup','Group']].head(2).reset_index(), 
            Reg.rename(columns={'index':'Coefficient'})], axis=1) ## Add drug name and group to regression output
        Out = Out.append(Reg)
    Out = Out.dropna()
    Out = Out[Out['Coefficient'] == 'LogScore']
    return(Out)


def GetPerDrugRegression():
    '''Performs linear regression on CCLE lines with load as the predictor and viability at the y.''' 
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset='CCLE', DataType='Drug', MutType='KsKa')
    Out = pd.DataFrame()
    for Drug in df['name'].unique():
        SingleDrug = ConvertPandasDFtoR(df[df['name'] == Drug].dropna())
        Reg= ConvertRDataframetoPandas(ro.r.DoLinearRegression(SingleDrug, NormalizeY=True)).reset_index()
        Reg = pd.concat([df.query('name == @Drug')[['name','subgroup','Group']].head(2).reset_index(), 
            Reg.rename(columns={'index':'Coefficient'})], axis=1) ## Add drug name and group to regression output
        Out = Out.append(Reg)
    return(Out)


def DoTTestForAS(row):
    '''This is a two-sided test for the null hypothesis that 2 independent samples have identical average (expected) values.'''
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
    Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
    Low = row[list(set(row.index) & set(Samples[Samples['MutLoad'] <= 10]['Barcode'].tolist()))].dropna()
    High = row[list(set(row.index) & set(Samples[Samples['MutLoad'] >= 1000]['Barcode'].tolist()))].dropna()
    return( pd.Series({'pVal': stats.ttest_ind(Low, High)[1], 
						'Low_mean': Low.mean(), 
						'High_mean': High.mean(), 
						'LowToHighDelta': Low.mean() - High.mean()}))


def GetDeltaPSIForEachGene():
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
    Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
    LowAndHigh = Samples[Samples['MutLoad'] < 10]['Barcode'].tolist() + Samples[Samples['MutLoad'] > 1000]['Barcode'].tolist()
    AS = GetAlternativeSplicingData(AS='RI', BarcodesOfInt=LowAndHigh)
    # Remove AS events if more than 25% of samples have a missing call
    PercentMissing = AS.isnull().sum(axis=1).reset_index()[0]/len(AS.columns)
    AS = AS.reset_index().assign(pct_missing = PercentMissing)
    AS = AS[AS['pct_missing'] < 0.25]
    # Get significantly different genes
    df = AS.apply(DoTTestForAS, axis=1).assign(GeneName=AS['symbol'])
    return(df)




# # GetNumberOfMutsAndCancerType(Dataset='TCGA')

# GetRegressionEstimates(Dataset='TCGA', DataType='AS', MutType='KsKa', GeneSet=[], ByGene=True).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ASMixedEffectRegressionEstimatesKsKaTCGA')


GetRegressionEstimates(Dataset='TCGA', DataType='Expression', MutType='KsKa', GeneSet=[], ByGene=True).to_csv(
    '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurityAndAge')

# GetRegressionEstimates(Dataset='CCLE', DataType='Expression', MutType='KsKa', GeneSet=[], ByGene=True, NormalizeY=True).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ExpressionOLSRegressionEstimatesKsKaCCLE')

# GetRegressionEstimates(Dataset='CCLE', DataType='Protein', MutType='KsKa', GeneSet=[], ByGene=True, NormalizeY=False).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ProteinOLSRegressionEstimatesKsKaCCLE')


# GetRegressionEstimates(Dataset='CCLE', DataType='RNAi', MutType='KsKa', GeneSet=[], ByGene=True, NormalizeY=False).to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/AchillesOnlyRNAiOLSRegressionEstimatesNoNormKsKaCCLE')


#RunVEPToGetDeleteriousScores(Dataset='CCLE').to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Raw/Mutations/CCLE/CCLE_Mutations_Polyphen_Annotated', sep='\t')



#GetHighlyExpressedGenes().to_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/HighlyExpressed_GTEX_GreaterThan500TPM.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='SNV').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_SNVs.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='Nonsynonymous').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_Nonsynonymous.txt')

# GetRegressionToGeneSetsOfInterest(Group='Chaperome', Dataset='TCGA', DataType='Expression', MutType='Polyphen').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/GroupedRegression_TCGA_Expression_Chaperome_Polyphen.txt')

