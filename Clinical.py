
import glob
import os
import pandas as pd
import numpy as np
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def SetUpPlottingPackages():
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


def GetCoxRegressionInput():
    '''
    Filtering patients for clinical analysis that:
    1) Have had no prior therapies (i.e. this treatment is the first line therapy)
    '''
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


def PlotCoxRegression(PlotType):
    SetUpPlottingPackages()
    df = GetCoxRegressionInput()
    if PlotType == 'ForestPlot':
        R = ConvertPandasDFtoR(df.astype(str))
        ro.r.GetForestPlot(R)
    elif PlotType == 'BroadLoad':
        df['MutBin'] = pd.cut(df['MutLoad'], bins=[0,50, 500,5600], labels=['0-50' ,'50-500','500-5000'])
        R = ConvertPandasDFtoR(df.astype(str))
        ro.r.PlotSurvivalCurve(R)



df = GetCoxRegressionInput()
df['MutBin'] = pd.cut(df['MutLoad'], bins=[0,50, 500,5600], labels=['0-50' ,'50-500','500-5000'])
R = ConvertPandasDFtoR(df.astype(str))



SetUpPlottingPackages()
ro.r.PlotSurvivalCurve(R)

SetUpPlottingPackages()
ro.r.GetForestPlot(R)

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

       