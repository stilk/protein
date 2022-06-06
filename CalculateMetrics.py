
import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import os.path
from GetRCodeIntoPython import *
from GetData import * 
from Splicing import * 


def GetRegressionStatsInput(Dataset, DataType, MutType, AllDrugs=False):
    ''' 
    Calculates mixed effects regression (controlling for cancer type) coefficients for
    the effect of expression/drug response/rnai response to mutational load.
    @Dataset = `CCLE`, `GTEX` or `TCGA`
    @DataType = `Expression`, `RNAi`, `Drug`, `AS` or `Protein`
    @AllDrugs = Optional variable for drug input data, if false will only report proteostasis drug
    '''
    if Dataset == "TS": # Just used for debugging; remove later
        return(GetExpressionData(Dataset='TS'))
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



def GetExpressionRegression(Dataset, DataType='Expression', MutType='KsKa', DoModelDiagnostics=False, ShuffleTMB=False):  
    SetUpRegressionPackages()
    df = GetRegressionStatsInput(Dataset, DataType, MutType)
    if ShuffleTMB: # Used to compare distribution of null vs empirical p-values for empirical FDR in Fig.1?
        df['LogScore'] = df.groupby(['type'])['LogScore'].transform(np.random.permutation) # shuffle mut load within cancer types 
    R = ConvertPandasDFtoR(df)
    if Dataset == 'CCLE':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'OLS', True, DoModelDiagnostics)))
    elif Dataset == 'TCGA':
        return(ConvertRDataframetoPandas(ro.r.DoRegressionPerGene(R, 'MixedEffect', True, DoModelDiagnostics)))
    

def DoRegressionByAge(DataType='Expression', MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/'):
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
    Output.to_csv(OutDir + 'TCGA_GLMM_ByAgeGroups')

def GetRegressionModelDiagnostics(Dataset, DataType, MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/Diagnostics/'):
    '''
    Calculates Pearson's R for fitted values with residuals for the full regression model, 
    and residuals vs all explanatory variables. Outputs statistics for each gene.
    '''
    if DataType == 'Expression':
        df = GetExpressionRegression(Dataset, DataType='Expression', MutType='KsKa', DoModelDiagnostics=True)
    elif DataType == 'shRNA':
        df = GetshRNARegression(DoModelDiagnostics=True)
    elif DataType == 'Drug':
        df = GetPerDrugRegression(DoModelDiagnostics=True)
    df.to_csv(OutDir + Dataset + DataType + MutType + 'ModelDiagnosticResidualVsFittedVals')
 


def JacknifeAcrossCancerTypes(Dataset, DataType, CancerTypeToRemove='', MutType='KsKa', OutDir=os.getcwd() + '/Data/Regression/Jacknife/'):
    ''' 
    Removes a cancer type from the dataset and re-calculates regression coefficients.
    Writes jacknifed regression coefficients to file since files per gene are huge.
    Only set up for per gene for now. 
    @CancerTypeToRemove = a string of which individual cancer type to remove. If empty loops through
        all cancer types. If not, removes cancer types of that string and exists loop after 1 iteration.
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
        


def GetshRNARegression(df = GetRegressionStatsInput(Dataset='CCLE', DataType='RNAi', MutType='KsKa'), DoModelDiagnostics=False):
    ''' Grouped regression '''
    SetUpRegressionPackages()
    groups = GetGeneAnnotationsOfInterest()
    df = df.merge(groups, left_on='GeneName', right_on='Hugo')
    Out = pd.DataFrame()
    for Group in groups['Group'].unique():
        SingleGroup = ConvertPandasDFtoR(df[df['Group'] == Group].dropna()[['LogScore','Value','GeneName']])
        #Reg= ConvertRDataframetoPandas(ro.r.DoLinearRegression(SingleGroup, NormalizeY=True)).reset_index()
        Reg = ConvertRDataframetoPandas(ro.r.DoOLSRegressionWithGeneName(SingleGroup, False, DoModelDiagnostics)).reset_index()
        Reg['Group'] = Group
        Reg['subgroup'] = groups[groups['Group'] == Group]['subgroup'].unique()[0]
        Out = Out.append(Reg)
    Out = Out.dropna().rename(columns={'index':'Coefficient'})
    if not DoModelDiagnostics:
        Out = Out[Out['Coefficient'] == 'LogScore']
    return(Out)



def GetPerDrugRegression(df = GetRegressionStatsInput(Dataset='CCLE', DataType='Drug', MutType='KsKa'), DoModelDiagnostics=False):
    '''Performs linear regression on CCLE lines with load as the predictor and viability at the y.''' 
    SetUpRegressionPackages()
    Out = pd.DataFrame()
    for Drug in df['name'].unique():
        SingleDrug = ConvertPandasDFtoR(df[df['name'] == Drug].dropna())
        Reg= ConvertRDataframetoPandas(ro.r.DoLinearRegression(SingleDrug, True, DoModelDiagnostics)).reset_index()
        Reg = pd.concat([df.query('name == @Drug')[['name','subgroup','Group','screenType']].head(2).reset_index(), 
            Reg.rename(columns={'index':'Coefficient'})], axis=1) ## Add drug name and group to regression output
        Out = Out.append(Reg)
    return(Out)


def DoTTestForAS(row):
    '''
    For each AS event, performs a two-sided t-test for the null hypothesis that the two groups 
    (low and high mutational load tumors) have the same expected average PSI values. 
    Outputs a Pandas Series of p-values for how different the two groups are and the average 
    change in PSI values between the two groups for each gene.
    Parameters 
    @row = a series of PSI values for each patient from @GetDeltaPSIForEachGene 
    '''
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
    Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
    Low = row[list(set(row.index) & set(Samples[Samples['MutLoad'] < 10]['Barcode'].tolist()))].dropna()
    High = row[list(set(row.index) & set(Samples[Samples['MutLoad'] > 1000]['Barcode'].tolist()))].dropna()
    return( pd.Series({'pVal': stats.ttest_ind(Low, High)[1], 
						'Low_mean': Low.mean(), 
						'High_mean': High.mean(), 
						'LowToHighDelta': Low.mean() - High.mean()}))


def GetDeltaPSIForEachGene():
    '''
    Performs a t-test for each intron retention event in TCGA between low and high mutational load tumors,
    removes intron retention events with missing calls and returns a table of average PSI changes for each gene,
    including if the change is significant.
    '''
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



def GetNumberGenesFilteredDueToPotentialEQTLs():
    '''
    Uses output from @GetExpLevelsForASEvents() to examine how many AS events are removed in
    both count categories (i.e., events splcied in/events not spliced in).
    '''
    ASDir='/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    Filt = pd.read_csv(ASDir +'TCGA_AS_EventsRemovedRI_Counts_ThresholdByPSI_0.8')
    Filt['Filtered'] = 1
    NotFilt= pd.read_csv(ASDir + 'TCGA_AS_EventsNOTRemovedRI_Counts_ThresholdByPSI_0.8')
    NotFilt['NotFiltered'] = 1
    FilteredStats = Filt.groupby(['Group'])['Filtered'].size().reset_index().merge(
        NotFilt.groupby(['Group'])['NotFiltered'].size().reset_index(),
        left_on=['Group'], right_on=['Group'])
    return(FilteredStats)




###########################################################
### Expression regression used for Fig. 2 TCGA and CCLE ###
###########################################################


# DoRegressionByAge()
# # GetNumberOfMutsAndCancerType(Dataset='TCGA')

# GetPerDrugRegression().to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/DrugOLSRegressionEstimatesKsKaCCLE')

# GetExpressionRegression(Dataset='TCGA', DataType='Expression', MutType='KsKa').to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ERChapAdded/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')

# GetCRISPRRegression().to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/CCLE_CRISPR_NormalizedByGene')


# GetshRNARegression().to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ERChapAdded/shRNA_GLMM_KsKa_CCLE')


# GetPerDrugRegression().to_csv(
#     '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/AllDrugsOLSRegressionEstimatesKsKaCCLE'
# # )

# GetExpressionRegression(Dataset='TCGA', DataType='Expression', MutType='KsKa', DoModelDiagnostics=False, ShuffleTMB=True).to_csv(
#    '/labs/ccurtis2/tilk/scripts/protein/Data/Regression/ERChapAdded/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurityShuffledLoad' 
# )

################
### Jacknife ###
################

# JacknifeAcrossCancerTypes('CCLE', DataType='RNAi', MutType='KsKa')

# JacknifeAcrossCancerTypes('CCLE', DataType='Drug', MutType='KsKa')

# JacknifeAcrossCancerTypes('TCGA', DataType='Expression', MutType='KsKa')

# JacknifeAcrossCancerTypes('TCGA', 'Expression', 'ACC')
