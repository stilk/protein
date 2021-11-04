
import glob
import os
import pandas as pd
import numpy as np
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def SetUpSurvivalPackages():
    ggplot = importr('ggplot2')
    #scales = importr('scales')
    #ggpmisc = importr('ggpmisc')
    cowplot = importr('cowplot')
    lme4 = importr('lme4')
    glmnet = importr('glmnet')
    lmertest = importr('lmerTest')
    data_table = importr('data.table')
    ro.r.source('/labs/ccurtis2/tilk/scripts/protein/Survival.R') # Source R script for plotting


def ConvertPandasDFtoR(df):
    with localconverter(ro.default_converter + pandas2ri.converter): 
        dfInR = ro.conversion.py2rpy(df) # Convert pandas dataframe to R 
    return(dfInR)


def ConvertRDataframetoPandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        dfInPd = ro.conversion.rpy2py(df) # Convert R dataframe back to Pandas
    return(dfInPd)

    pd.concat([
        pd.DataFrame({'Group': '20S Core' , 'subgroup' : 'Proteasome', 'Hugo' : corum[corum['ComplexName'] == '20S proteasome']['subunits(Gene name)'].str.split(';', expand=True).melt()['value']  }),
        pd.DataFrame({'Group': '19S Regulatory Particle' , 'subgroup' : 'Proteasome', 'Hugo' : ['PSMC1','PSMC2','PSMC3','PSMC4','PSMC5','PSMC6', 'PSMD1','PSMD2','PSMD3','PSMD4','PSMD6','PSMD7','PSMD8','PSMD11','PSMD12','PSMD13','PSMD14'] })
    ])

def GetTCGAClinical():
    DataDir = os.getcwd() + '/Data/Raw/Clinical/'
    # Get clinical info and match to TMB
    clin = pd.read_csv(DataDir + '/TCGA/TCGA-CDR-SupplementalTableS1.txt', sep='\t')
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
    Samples['Barcode'] = Samples['Barcode'].str[0:12]
    clin = clin.merge(Samples, left_on='bcr_patient_barcode', right_on='Barcode')
    # Get average gene expression of the proteasome and match for all cancer types 
    exp = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_AvgGeneExpForAllCancerTypes')
    exp['Barcode'] = exp['Barcode'].str[0:12]
    corum = pd.read_csv('/labs/ccurtis2/tilk/02_DNDS/separateDatabaseFiles/CORUM/coreComplexes.txt.zip', sep='\t')
    Proteasome = pd.concat([
        pd.DataFrame({'Group': '20S Core' , 'subgroup' : 'Proteasome', 'Hugo' : corum[corum['ComplexName'] == '20S proteasome']['subunits(Gene name)'].str.split(';', expand=True).melt()['value']  }),
        pd.DataFrame({'Group': '19S Regulatory Particle' , 'subgroup' : 'Proteasome', 'Hugo' : ['PSMC1','PSMC2','PSMC3','PSMC4','PSMC5','PSMC6', 'PSMD1','PSMD2','PSMD3','PSMD4','PSMD6','PSMD7','PSMD8','PSMD11','PSMD12','PSMD13','PSMD14'] })
    ])
    # Set a threshold that if half of proteasomal genes are over-expressed, we'll call this 'over-expressed' and vice versa
    exp=exp.merge(Proteasome, left_on='GeneName', right_on='Hugo')
    exp.to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_AvgGeneExpInProteasomeForAllCancerTypes')
    thresh = exp.groupby(['Barcode','Group'])[['DepletedTranscript','OverExpressedTranscript']].sum().reset_index()
    thresh['LowProteasome'] = thresh['DepletedTranscript'] > len(Proteasome)/2
    thresh['HighProteasome'] = thresh['OverExpressedTranscript'] > len(Proteasome)/2
    clin = clin.merge(thresh, left_on='bcr_patient_barcode', right_on='Barcode')
    clin.to_csv(DataDir + 'AnnotatedTCGAClinicalOutcomes.txt')


def GetMMRFData(DataType):
    '''
    @DataType = string of either `Mutations` or `Expression`
    '''
    ListOfFiles = glob.glob(os.path.join(os.getcwd() + '/Data/Raw/' + DataType + '/MMRF/*/', "*gz"))
    Out = pd.DataFrame()
    for FileName in ListOfFiles:
        if DataType == 'Expression':
            Out = Out.append(pd.read_csv(FileName, sep='\t', header=None))
        elif DataType == 'Mutations':
            Muts = pd.read_csv(FileName, skiprows=7, sep='\t')
            Muts = Muts[~Muts['Variant_Classification'].isin(['In_Frame_Ins','In_Frame_Del'])] # Only look at SNVs
            Muts = Muts.groupby(['Tumor_Sample_Barcode']).size().reset_index().rename(columns={0:'MutLoad'})
            Out = Out.append(Muts)
        elif DataType == 'Clinical':
            Out = Out.append(pd.read_csv(FileName, sep='\t'))
    return(Out)


def GetCoxRegressionInput(Dataset):
    '''
    Filtering patients for clinical analysis that:
    1) Have had no prior therapies (i.e. this treatment is the first line therapy)
    '''
    if Dataset == 'MMRF':
        muts = GetMMRFData('Mutations')
        muts['Barcode'] = muts['Tumor_Sample_Barcode'].str[0:9]
        df = GetMMRFData(DataType='Clinical')
        df = df.query("regimen_or_line_of_therapy == 'First line of therapy'")
        df = df.replace({"'--" : np.nan}).dropna(axis=1, how='all') # remove columns that don't have info
        df = df.merge(muts, left_on='case_submitter_id', right_on='Barcode')
        NumTreatPerPatient = df.groupby(['case_submitter_id'])['therapeutic_agents'].size().reset_index().rename(
            columns={'therapeutic_agents':'NumberOfDiffTreatments'})
        MetaData = df[['iss_stage','gender','case_submitter_id','age_at_index','days_to_last_follow_up','vital_status']].drop_duplicates().rename(
            columns={'gender':'sex'})
        MetaData['survival'] = MetaData['vital_status'].replace('Alive',0).replace('Dead',1) # recode survival as 0 or 1
        MetaData = MetaData.merge(NumTreatPerPatient, left_on='case_submitter_id', right_on='case_submitter_id')
        df['ID'] = 1
        pivoted = pd.pivot_table(df, values='ID', index=['case_submitter_id','MutLoad'], columns = 'therapeutic_agents').reset_index().replace(np.nan,0)
        # Has to have received one of the 3 proteasome inhibitors in this cohort
        pivoted = pivoted[(pivoted['Bortezomib'] > 0) | (pivoted['Carfilzomib'] > 0) | (pivoted['Ixazomib'] > 0)]
        pivoted = pivoted.merge(MetaData, left_on='case_submitter_id', right_on='case_submitter_id')
        return(pivoted)
    elif Dataset == 'TCGA':
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein//Data/Raw/Clinical/AnnotatedTCGAClinicalOutcomes.txt')
        df = df[df['Group'] == '19S Regulatory Particle']
        df = df[df['Redaction'] != 'Redacted']
        df['HighProteasome'] = (df['OverExpressedTranscript'] > 0) & (df['DepletedTranscript'] == 0)
        df['LowProteasome'] = (df['OverExpressedTranscript'] == 0) & (df['DepletedTranscript'] > 0)
        #df['OverExpressedFraction'] = df['OverExpressedTranscript']/ (df['OverExpressedTranscript'] + df['DepletedTranscript'])
        #df['DepletedFraction'] = df['DepletedTranscript']/ (df['OverExpressedTranscript'] + df['DepletedTranscript'])
        #df['HighProteasome'] = df['OverExpressedFraction'] > 0.5
        #df['LowProteasome'] = df['DepletedFraction'] > 0.5
        # df['HighProteasome'] = (df['OverExpressedTranscript'] > 0))
        # df['LowProteasome'] = ( df['DepletedTranscript'] > 0 ) & (df['OverExpressedTranscript'] == 0)
        df['MutLoadGroup'] = np.where(df['MutLoad'] <= 10, 'Low', np.NaN)
        df['MutLoadGroup'] = np.where(df['MutLoad'] >= 1000 , 'High', df['MutLoadGroup'])
        #df = df[(df['LowProteasome'] == True) | (df['HighProteasome'] == True)]
        df['ProteasomeGroup'] = np.where((df['LowProteasome'] == True), 'LowProteasome', '')
        df['ProteasomeGroup'] = np.where((df['HighProteasome'] == True), 'HighProteasome', df['ProteasomeGroup'])
        df['StrataGroup'] = np.where((df['LowProteasome'] == True) & (df['MutLoadGroup'] == 'Low'), 'LowPro_LowTMB','')
        df['StrataGroup'] = np.where((df['LowProteasome'] == True) & (df['MutLoadGroup'] == 'High'), 'LowPro_HighTMB', df['StrataGroup'])
        df['StrataGroup'] = np.where((df['HighProteasome'] == True) & (df['MutLoadGroup'] == 'Low'), 'HighPro_LowTMB',df['StrataGroup'])
        df['StrataGroup'] = np.where((df['HighProteasome'] == True) & (df['MutLoadGroup'] == 'High'), 'HighPro_HighTMB',df['StrataGroup'])
        df = df[~(df['StrataGroup'] == '')]
        return(df)


test = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_AvgGeneExpInProteasomeForAllCancerTypes')
test2 = pd.concat([test.groupby(['Barcode'])['DepletedTranscript'].value_counts().unstack().reset_index().replace(np.nan,0).assign(group='DepletedTranscript'),
                   test.groupby(['Barcode'])['OverExpressedTranscript'].value_counts().unstack().reset_index().replace(np.nan,0).assign(group='OverExpressedTranscript')])

Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
Samples['Barcode'] = Samples['Barcode'].str[0:12]
test3 = test2.merge(Samples, left_on='Barcode', right_on='Barcode')


def PlotCoxRegression(PlotType, Dataset):
    SetUpSurvivalPackages()
    df = GetCoxRegressionInput(Dataset)
    if PlotType == 'ForestPlot':
        R = ConvertPandasDFtoR(df.astype(str))
        SetUpSurvivalPackages();ro.r.GetForestPlot(R, Dataset)
        ro.r.GetForestPlot(R, Dataset)
    elif PlotType == 'BroadLoad':
        df['MutBin'] = pd.cut(df['MutLoad'], bins=[0,50, 500,5600], labels=['0-50' ,'50-500','500-5000'])
        R = ConvertPandasDFtoR(df.astype(str))
        ro.r.PlotSurvivalCurve(R, Dataset)




GetTCGAClinical()

# R =ConvertPandasDFtoR(GetCoxRegressionInput(Dataset).astype(str))
# SetUpSurvivalPackages();ro.r.PlotSurvivalCurve(R, Dataset)


# R =ConvertPandasDFtoR(GetCoxRegressionInput(Dataset).astype(str))
# SetUpSurvivalPackages();ro.r.GetForestPlot(R, Dataset)


# m=GetCoxRegressionInput(Dataset='MMRF')
# test= m[(m['NumberOfDiffTreatments'] <= 3)]
# test['MutLoadGroup'] = np.where(test['MutLoad'] < 100, 'Low', np.NaN)
# test['MutLoadGroup'] = np.where(test['MutLoad'] >= 100 , 'High', test['MutLoadGroup'])
# R = ConvertPandasDFtoR(test.astype(str))
# SetUpSurvivalPackages();ro.r.PlotSurvivalCurve(R, 'MMRF')

# # Remove drugs that have opposite effect of PI that induce 
# # pivoted = pivoted[(pivoted['Lenalidomide'] == 0) & (pivoted['Thalidomide'] == 0) & (pivoted['Pomalidomide'] == 0)]


# PI = ['Bortezomib', 'Carfilzomib','Ixazomib']

#  lenalidomide and pomalidomide
# ', '','P' ## INDUCES DEGRADATION BY DRUG

# Immuno = ['Elotuzumab','Daratumumab','Dexamethasone','Pomalidomide','Prednisone','Thalidomide']
# Chemo = ['Bendamustine','Cyclophosphamide', 'Doxorubicin','Lenalidomide','Melphalan']

    

# pivoted = pd.pivot_table(df, values='ID', 
#  index=['iss_stage','gender','race','case_submitter_id','age_at_index','days_to_death','vital_status'], 
#  columns = 'therapeutic_agents').reset_index()

# pivoted = pd.pivot_table(df, values='ID', 
#  index=['case_submitter_id'], 
#  columns = 'therapeutic_agents').reset_index().np.nan(0)


# 'case_submitter_id', 'days_to_treatment_end','days_to_treatment_start','therapeutic_agents','days_to_death', 'gender', 'age_at_diagnosis'

       