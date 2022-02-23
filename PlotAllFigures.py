
import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import FloatVector
from CalculateMetrics import *


def SetUpPlottingPackages():
    ggplot = importr('ggplot2')
    #scales = importr('scales')
    #ggpmisc = importr('ggpmisc')
    cowplot = importr('cowplot')
    lme4 = importr('lme4')
    glmnet = importr('glmnet')
    lmertest = importr('lmerTest')
    data_table = importr('data.table')
    dplyr = importr('dplyr')
    gprofiler = importr('gprofiler2')
    ro.r.source('/labs/ccurtis2/tilk/scripts/protein/Plotting.R') # Source R script for plotting
    ro.r.source('/labs/ccurtis2/tilk/scripts/protein/DoGeneSetEnrichment.R') # Source gse script



def ConvertPandasDFtoR(df):
    with localconverter(ro.default_converter + pandas2ri.converter): 
        dfInR = ro.conversion.py2rpy(df) # Convert pandas dataframe to R 
    return(dfInR)


def ConvertRDataframetoPandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        dfInPd = ro.conversion.rpy2py(df) # Convert R dataframe back to Pandas
    return(dfInPd)


def GetGeneAnnotationsOfInterest():
    corum = corum = pd.read_csv('/labs/ccurtis2/tilk/02_DNDS/separateDatabaseFiles/CORUM/coreComplexes.txt.zip', sep='\t')
    chap = pd.read_csv('/labs/ccurtis2/tilk/09_PROTEIN/Human_Chaperome_TableS1A_PMID_29293508', sep='\t')
    groups = pd.concat([
        pd.DataFrame({'Group': 'Cytoplasmic Ribosomes','subgroup' : 'Ribosomes', 'Hugo' : corum[corum['ComplexName'].str.contains('subunit, cytoplasmic')]['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].drop_duplicates().tolist()}),
        pd.DataFrame({'Group': 'Mitochondrial Ribosomes', 'subgroup' : 'Ribosomes', 'Hugo' : corum[corum['ComplexName'].str.contains('subunit, mitochondrial')]['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].drop_duplicates().tolist()}),
        pd.DataFrame({'Group': '20S Core' , 'subgroup' : 'Proteasome', 'Hugo' : corum[corum['ComplexName'] == '20S proteasome']['subunits(Gene name)'].str.split(';', expand=True).melt()['value']  }),
        pd.DataFrame({'Group': '19S Regulatory Particle' , 'subgroup' : 'Proteasome', 'Hugo' : ['PSMC1','PSMC2','PSMC3','PSMC4','PSMC5','PSMC6', 'PSMD1','PSMD2','PSMD3','PSMD4','PSMD6','PSMD7','PSMD8','PSMD11','PSMD12','PSMD13','PSMD14'] }),
        pd.DataFrame({'Group': 'Mitochondrial Chaperones' , 'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'MITO') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'ER Chaperones' , 'subgroup' : 'Chaperones','Hugo' : chap[(chap['Level2'] == 'ER') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'Small HS' , 'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'sHSP') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 100' ,'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'HSP100') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 90', 'subgroup' : 'Chaperones', 'Hugo' :chap[(chap['Level2'] == 'HSP90') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 70','subgroup' : 'Chaperones', 'Hugo' : chap[(chap['Level2'] == 'HSP70') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 60','subgroup' : 'Chaperones', 'Hugo' : chap[(chap['Level2'] == 'HSP60') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']}),
        pd.DataFrame({'Group': 'HSP 40', 'subgroup' : 'Chaperones','Hugo' : chap[(chap['Level2'] == 'HSP40') & (chap['CHAP / CO-CHAP'] == 'CHAP') ]['Gene']})])
    return(groups)


def GetFigureInput(FigureNum):
    SetUpPlottingPackages()
    rstats = importr('stats')
    InputDir = '/labs/ccurtis2/tilk/scripts/protein/Data/'
    if FigureNum == 'GlobalGSE_TCGA_Regression': # fig1B and 1C
        df = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['GeneName'] = df['GeneName'].str.split('_', expand=True)[0]
        df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['Pr...t..']), method = 'fdr')
        GeneSet = df[(df['Estimate'] > 0) & (df['Adj.Pval'] < 0.05)]['GeneName']
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(GeneSet)))
        return(ConvertPandasDFtoR(gse[['term_name','source','p_value']]))
    elif FigureNum == 'Groups_CCLEAndTCGA':
        tcga = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGA').assign(Dataset='TCGA')
        ccle = pd.read_csv(InputDir + 'Regression/ExpressionOLSRegressionEstimatesKsKaCCLE').assign(Dataset='CCLE')
        all = pd.concat([tcga[['Pr...t..','adj.pval','Estimate','GeneName','Unnamed: 0','Dataset']], 
                         ccle[['Pr...t..','adj.pval','Estimate','GeneName','Unnamed: 0','Dataset']]]).rename(columns={'Pr...t..':'pval'})
        all['Coefficient'] = all['Unnamed: 0'].str.replace('\d+', '')
        all = all[all['Coefficient'] == 'LogScore']
        quantile = all.groupby(['Dataset'])['Estimate'].quantile([0,0.1,0.5,0.9,1]).reset_index().rename(columns={'level_1':'Group'}).assign(
        subgroup='Quantile').assign(pval = 0).assign(GeneName='')
        groups = GetGeneAnnotationsOfInterest()
        all = all.merge(groups, left_on='GeneName', right_on='Hugo')
        all = pd.concat([all[['Dataset','Group','Estimate','subgroup','pval','GeneName']], quantile])
        return(ConvertPandasDFtoR(all.astype(str))) 
    elif FigureNum == 'AS_Delta_PSI': # Fig 2
        SetUpPlottingPackages()
        df = GetDeltaPSIForEachGene()
        df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['pVal']), method = 'fdr')
        NegGenes = df[(df['Adj.Pval'] < 0.05) & (df['LowToHighDelta'] < -0.01)]['GeneName'].drop_duplicates()
        PosGenes = df[(df['Adj.Pval'] < 0.05) & (df['LowToHighDelta'] > 0.01)]['GeneName'].drop_duplicates()
        NegGSE = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(NegGenes)))
        PosGSE = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(PosGenes)))
        out = pd.concat([PosGSE.assign(Group= 'PosPSI'), NegGSE.assign(Group= 'NegPSI')])
        out['ID'] = np.arange(0, len(out))
        out = out[out['source'].isin(['CORUM','REAC'])].set_index(['term_name']).groupby(
            ['Group','source'])['p_value'].nsmallest(10).reset_index()
        return(ConvertPandasDFtoR(out))
    elif FigureNum == 'JacknifeExpCCLE':
        Complexes = GetSubsetOfGeneAnnotationsOfInterest()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLEExpressionKsKa*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df[df['Coefficient'] == 'LogScore'] 
            df = df.merge(Complexes, right_on='Hugo', left_on='GeneName')
            print(df)
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureNum == 'JacknifeExpTCGA':
        Complexes = GetSubsetOfGeneAnnotationsOfInterest()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/TCGAExpression*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df[df['Coefficient'] == 'LogScore'] 
            df = df.merge(Complexes, right_on='Hugo', left_on='GeneName')
            print(df)
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureNum == 'JacknifeDrugsCCLE': #fig 3
        DrugsOfInterest = GetDrugResponseData(Screen='primary', AllDrugs=False)[['name','Group','subgroup']].drop_duplicates()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLEDrug*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df = df[df['Coefficient'] == 'LogScore'] 
            Out = Out.append(df.drop(columns={'Unnamed: 0'}))
        df = GetPerDrugRegression().assign(CancerTypeRemoved='All Cancers')
        Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureNum == 'JacknifeshRNACCLE':
        Complexes = GetSubsetOfGeneAnnotationsOfInterest()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLERNAi*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df = df[df['Coefficient'] == 'LogScore'] 
            Out = Out.append(df.drop(columns={'Unnamed: 0'}))
        df = GetshRNARegression().assign(CancerTypeRemoved='All Cancers')
        Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureNum == 'AllDrugsByLoad':
        df = pd.read_csv(os.getcwd() + '/Data/Regression/AllDrugsOLSRegressionEstimatesKsKaCCLE')
        df = df[df['Coefficient'] == 'LogScore'] 
        return(ConvertPandasDFtoR(df))
    elif FigureNum == 'Multicollinearity':
        from functools import reduce
        ccle_muts = AnnotateMutationalLoad(GetPointMutations('CCLE'), 'SNV').rename(columns={'MutLoad':'SNV'}).merge(
                AnnotateMutationalLoad(GetPointMutations('CCLE'), 'KsKa').rename(columns={'MutLoad':'Protein-Coding'}),
                left_on='Barcode',right_on='Barcode',how='outer')
        ccle_cnvs = AnnotateMutationalLoad(GetCopyNumberMutations('CCLE'), 'CN_Deletions').rename(columns={'MutLoad':'Deletions'}).merge(
                AnnotateMutationalLoad(GetCopyNumberMutations('CCLE'), 'CN_Amplifications').rename(columns={'MutLoad':'Amplifications'}), 
                left_on='Barcode',right_on='Barcode',how='outer')
        ccle_merged = ccle_muts.merge(ccle_cnvs, left_on='Barcode',right_on='Barcode',how='outer').dropna()
        tcga_muts = AnnotateMutationalLoad(GetPointMutations('TCGA'), 'SNV').rename(columns={'MutLoad':'SNV'}).merge(
                AnnotateMutationalLoad(GetPointMutations('TCGA'), 'KsKa').rename(columns={'MutLoad':'Protein-Coding'}), 
                left_on='Barcode',right_on='Barcode',how='outer')
        tcga_cnvs = AnnotateMutationalLoad(GetCopyNumberMutations('TCGA'), 'CN_Deletions').rename(columns={'MutLoad':'Deletions'}).merge(
                AnnotateMutationalLoad(GetCopyNumberMutations('TCGA'), 'CN_Amplifications').rename(columns={'MutLoad':'Amplifications'}),
                left_on='Barcode',right_on='Barcode',how='outer')
        tcga_merged = tcga_muts.merge(tcga_cnvs, left_on='Barcode',right_on='Barcode',how='outer').dropna()  
        return(ConvertPandasDFtoR(pd.concat([tcga_merged.assign(Dataset='TCGA'), ccle_merged.assign(Dataset='CCLE')])))



def GetFigure(Figure):
    if Figure == 'Groups_CCLEAndTCGA': # Fig 2A
        all = GetFigureInput('Groups_CCLEAndTCGA')
        SetUpPlottingPackages(); ro.r.PlotRegCoefPerGroup(all)
    elif Figure == 'AS_Delta_PSI': # Fig 2B
        out = GetFigureInput('AS_Delta_PSI')
        SetUpPlottingPackages(); ro.r.PlotDeltaPSI(out)
    elif Figure == 'AS_PSI': #Fig 2C
        SetUpPlottingPackages(); ro.r.VisualizeAS()
    elif Figure == 'GlobalGSE_TCGA_Regression':  #fig 1B
        foo = GetFigureInput('GlobalGSE_TCGA_Regression')
        SetUpPlottingPackages(); ro.r.PlotCircularCORUMTCGA(foo)
    elif Figure == 'GlobalGSE_TCGA_Regression_KEGG': #fig 1C
        foo = GetFigureInput('GlobalGSE_TCGA_Regression')
        SetUpPlottingPackages(); ro.r.PlotCircularKeggTCGA(foo)
    elif Figure == 'JacknifeExpCCLE': # supplemental figure x
        df = GetFigureInput('JacknifeExpCCLE')
        SetUpPlottingPackages(); ro.r.PlotJacknifedExpressionAcrossGroups(df, 'CCLE')
    elif Figure == 'JacknifeExpTCGA':  # supplemental figure x
        df = GetFigureInput('JacknifeExpTCGA')
        SetUpPlottingPackages(); ro.r.PlotJacknifedExpressionAcrossGroups(df,'TCGA')
    elif Figure == 'JacknifeDrugsCCLE': # fig 3b
        df = GetFigureInput('JacknifeDrugsCCLE')
        SetUpPlottingPackages(); ro.r.PlotJacknifedDrugsAcrossGroups(df)
    elif Figure == 'JacknifeshRNACCLE': # fig 3a
        df = GetFigureInput('JacknifeshRNACCLE')
        SetUpPlottingPackages(); ro.r.PlotJacknifedshRNAAcrossGroups(df)
    elif Figure == 'AllDrugsByLoad': # fig4a
        df = GetFigureInput('AllDrugsByLoad')
        SetUpPlottingPackages(); ro.r.PlotFractionOfSigDrugsForAll(df)
        #SetUpPlottingPackages(); ro.r.PlotRegressionForAllDrugs(df)
    elif Figure == 'BootstrappedDrugs': # fig4b
        SetUpPlottingPackages(); ro.r.PlotBootstrappedNegativelyAssociatedDrugs()
    elif Figure == 'Multicollinearity': # supplemental figure x
        df = GetFigureInput('Multicollinearity')
        SetUpPlottingPackages(); ro.r.PlotMultiCollinearityOfSNVsAndCNVs(df)
    elif Figure == 'PSI_Supplemental':
        ro.r.VisualizeAllASThresholds() # supplemental figX A
        ro.r.VisualizeAllAS() # supplemental figX B
        

    



    # elif Figure == 'DeltaPSI':     
    #     df['pVal'] < 0.05]
    #     gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ro.r.GetGeneSets('ExpressionCCLE')))
    #     gse = gse[gse['source'] == 'CORUM']
    #     return(ConvertPandasDFtoR(gse[['term_name','p_value']]))




### GENE SET ENRICHMENT PLOTS -- FIG 1
# SetUpPlottingPackages(); ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_TCGA_Regression'))
# SetUpPlottingPackages(); ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_CCLE_Regression'))
# SetUpPlottingPackages(); ro.r.PlotGlobalGSETCGA(GetFigureInput('DGE_GSE_TCGA'))

### EXPLORING BY CATEGORY -- FIG 2
# SetUpPlottingPackages(); GetFigure(Figure='Groups_TCGA')
# SetUpPlottingPackages(); GetFigure(Figure='Groups_CCLE')
# SetUpPlottingPackages(); GetFigure(Figure='Protein_Exp')


GetGLMMDiagnostics()