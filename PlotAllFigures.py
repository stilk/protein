'''
This script runs all of the analyses and plots all of the figures in the manuscript.
'''

import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
from CalculateMetrics import *
from Splicing import *
from GetAnnotations import *


def GetFigureInput(FigureName):
    '''
    Outputs a dataframe in R of all of the raw data used for plotting each supplemental and main figure in the manuscript.
    Figures are all plotted using @GetFigure. Commented out functions 
    Parameters:
        @FigureName = a unique string identifier for each figure
    '''
    rstats = importr('stats')
    Complexes = GetGeneAnnotationsOfInterest()
    InputDir = os.getcwd() + '/Data/'
    if FigureName == 'GlobalGSE_TCGA_Regression': # Data for Fig 1B and 1C
        #GetExpressionRegression(Dataset='TCGA', DataType='Expression', MutType='KsKa').to_csv(
        #   InputDir + 'ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')
        df = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['GeneName'] = df['GeneName'].str.split('_', expand=True)[0]
        df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['Pr...t..']), method = 'fdr')
        GeneSet = df[(df['Estimate'] > 0) & (df['Adj.Pval'] < 0.05)]['GeneName']
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(GeneSet)))
        return(ConvertPandasDFtoR(gse[['term_name','source','p_value']]))
    elif FigureName == 'AS_PSI': # Data for 2A
        # GetExpLevelsForASEvents(ASType='RI', PSI_Threshold=0.8, FilterForeQTLS=True)
        df = pd.read_csv(InputDir + '/AS_Tables/TCGA_RI_Counts_ThresholdByPSI_0.8')
        df['Filtered'] = 'Filtered'
        return(ConvertPandasDFtoR(df))
    elif FigureName == 'AS_Delta_PSI': # Data for 2B
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
    elif FigureName == 'Groups_CCLEAndTCGA': # Data for Fig 3
        #GetExpressionRegression(Dataset='TCGA', DataType='Expression', MutType='KsKa').to_csv(
        #   InputDir + 'ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')
        tcga = pd.read_csv(InputDir + 'Regression/ERChapAdded/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity').assign(Dataset='TCGA')
        #GetExpressionRegression(Dataset='CCLE', DataType='Expression', MutType='KsKa').to_csv(
        #   InputDir + 'ExpressionOLSRegressionEstimatesKsKaCCLE')
        ccle = pd.read_csv(InputDir + 'Regression/ERChapAdded/ExpressionOLSRegressionEstimatesKsKaCCLE').assign(Dataset='CCLE')
        all = pd.concat([tcga[['Pr...t..','Estimate','GeneName','Unnamed: 0','Dataset']], 
                         ccle[['Pr...t..','Estimate','GeneName','Unnamed: 0','Dataset']]]).rename(columns={'Pr...t..':'pval'})
        all['Coefficient'] = all['Unnamed: 0'].str.replace('\d+', '')
        all = all[all['Coefficient'] == 'LogScore']
        quantile = all.groupby(['Dataset'])['Estimate'].quantile([0,0.1,0.5,0.9,1]).reset_index().rename(columns={'level_1':'Group'}).assign(
            subgroup='Quantile').assign(pval = 0).assign(GeneName='')
        groups = GetGeneAnnotationsOfInterest()
        all = all.merge(groups, left_on='GeneName', right_on='Hugo')
        all = pd.concat([all[['Dataset','Group','Estimate','subgroup','pval','GeneName']], quantile])
        return(ConvertPandasDFtoR(all.astype(str))) 
    elif FigureName == 'JacknifeshRNACCLE': # Data for Fig 4A
        # JacknifeAcrossCancerTypes('CCLE', DataType='RNAi', MutType='KsKa')
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLERNAi*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df = df[df['Coefficient'] == 'LogScore'] 
            Out = Out.append(df.drop(columns={'Unnamed: 0'}))
        df = GetshRNARegression().assign(CancerTypeRemoved='All Cancers')
        Out = Out[['Estimate','Pr...t..','Group','subgroup','CancerTypeRemoved']].append(df[['Estimate','Pr...t..','Group','subgroup','CancerTypeRemoved']])
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'JacknifeDrugsCCLE': # Data for Fig 4B
        # JacknifeAcrossCancerTypes('CCLE', DataType='Drug', MutType='KsKa')
        DrugsOfInterest = GetDrugResponseData(Screen='primary', AllDrugs=False)[['name','Group','subgroup']].drop_duplicates()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLEDrug*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df = df[df['Coefficient'] == 'LogScore'] 
            Out = Out.append(df.drop(columns={'Unnamed: 0'}))
        df = GetPerDrugRegression().assign(CancerTypeRemoved='All Cancers')
        df = df[df['Coefficient'] == 'LogScore'] 
        Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'AllDrugsByLoad': # Data for Fig 5A
        # GetPerDrugRegression()
        df = pd.read_csv(os.getcwd() + '/Data/Regression/AllDrugsOLSRegressionEstimatesKsKaCCLE')
        df = df[df['Coefficient'] == 'LogScore'] 
        return(ConvertPandasDFtoR(df))
    elif FigureName == 'Multicollinearity': # Data for Sup. Fig 1
        ccle_muts = AnnotateMutationalLoad(GetPointMutations('CCLE'), 'SNV').rename(columns={'MutLoad':'SNV'}).merge(
                    AnnotateMutationalLoad(GetPointMutations('CCLE'), 'KsKa').rename(columns={'MutLoad':'Protein-Coding'}),
                    left_on='Barcode',right_on='Barcode',how='outer').merge(AnnotateMutationalLoad(GetPointMutations('CCLE'), 'Synonymous'),
                    left_on='Barcode',right_on='Barcode',how='outer').rename(columns={'MutLoad':'Synonymous'}).merge(
                    AnnotateMutationalLoad(GetPointMutations('CCLE'), 'Nonsynonymous'),left_on='Barcode',right_on='Barcode' ,how='outer').rename(
                    columns={'MutLoad':'Nonsynonymous'}).merge( AnnotateMutationalLoad(GetPointMutations('CCLE'), 'Nonsense'),
                    left_on='Barcode',right_on='Barcode',how='outer').rename(columns={'MutLoad':'Nonsense'})
        ccle_cnvs = AnnotateMutationalLoad(GetCopyNumberMutations('CCLE'), 'CN_Deletions').rename(columns={'MutLoad':'Deletions'}).merge(
                    AnnotateMutationalLoad(GetCopyNumberMutations('CCLE'), 'CN_Amplifications').rename(columns={'MutLoad':'Amplifications'}), 
                    left_on='Barcode',right_on='Barcode',how='outer').merge(AnnotateMutationalLoad(GetCopyNumberMutations('CCLE'), 'CNV'),
                    left_on='Barcode',right_on='Barcode',how='outer').rename(columns={'MutLoad':'CNV'})
        ccle_merged = ccle_muts.merge(ccle_cnvs, left_on='Barcode',right_on='Barcode',how='outer').dropna()
        tcga_muts = AnnotateMutationalLoad(GetPointMutations('TCGA'), 'SNV').rename(columns={'MutLoad':'SNV'}).merge(
                    AnnotateMutationalLoad(GetPointMutations('TCGA'), 'KsKa').rename(columns={'MutLoad':'Protein-Coding'}),
                    left_on='Barcode',right_on='Barcode',how='outer').merge(AnnotateMutationalLoad(GetPointMutations('TCGA'), 'Synonymous'),
                    left_on='Barcode',right_on='Barcode',how='outer').rename(columns={'MutLoad':'Synonymous'}).merge(
                    AnnotateMutationalLoad(GetPointMutations('TCGA'), 'Nonsynonymous'),left_on='Barcode',right_on='Barcode' ,how='outer').rename(
                    columns={'MutLoad':'Nonsynonymous'}).merge( AnnotateMutationalLoad(GetPointMutations('TCGA'), 'Nonsense'),
                    left_on='Barcode',right_on='Barcode',how='outer').rename(columns={'MutLoad':'Nonsense'})
        tcga_cnvs = AnnotateMutationalLoad(GetCopyNumberMutations('TCGA'), 'CN_Deletions').rename(columns={'MutLoad':'Deletions'}).merge(
                    AnnotateMutationalLoad(GetCopyNumberMutations('TCGA'), 'CN_Amplifications').rename(columns={'MutLoad':'Amplifications'}), 
                    left_on='Barcode',right_on='Barcode',how='outer').merge(AnnotateMutationalLoad(GetCopyNumberMutations('TCGA'), 'CNV'),
                    left_on='Barcode',right_on='Barcode',how='outer').rename(columns={'MutLoad':'CNV'})
        tcga_merged = tcga_muts.merge(tcga_cnvs, left_on='Barcode',right_on='Barcode',how='outer').dropna()  
        return(ConvertPandasDFtoR(pd.concat([tcga_merged.assign(Dataset='TCGA'), ccle_merged.assign(Dataset='CCLE')])))
    elif FigureName == 'NumASEventsFiltered':  # Data for Sup. Fig 2A
        filtered = GetNumberGenesFilteredDueToPotentialEQTLs()
        return(ConvertPandasDFtoR(filtered))
    elif FigureName == 'AS_PSI_FilteredCounts': # Data for Sup. Fig 2B
        # GetExpLevelsForASEvents(ASType='RI', PSI_Threshold=0.8, FilterForeQTLS=False)
        df = pd.read_csv(InputDir + '/AS_Tables/TCGA_RI_Counts_ThresholdByPSI_0.8FiltereQTLS_False')
        df['Filtered'] = 'NotFiltered'
        return(ConvertPandasDFtoR(df))
    elif FigureName == 'JacknifeExpTCGA': # Data for Sup. Fig 4
        # JacknifeAcrossCancerTypes('TCGA', DataType='Expression', MutType='KsKa')
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/TCGAExpression*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df[df['Coefficient'] == 'LogScore'] 
            df = df.merge(Complexes, right_on='Hugo', left_on='GeneName')
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'CorrelationWithAge': # Data for Sup. Fig 5
        # DoRegressionByAge()
        exp = pd.read_csv(InputDir + '/Regression/TCGA_GLMM_ByAgeGroups')
        exp['Coefficient'] = exp['Unnamed: 0'].str.replace('\d+', '')
        exp = exp[exp['Coefficient'] == 'LogScore']
        exp = exp.merge(GetGeneAnnotationsOfInterest(), left_on='GeneName', right_on='Hugo', how='left')
        return(ConvertPandasDFtoR(exp.dropna())) 
    elif FigureName == 'JacknifeExpCCLE': # Data for Sup. Fig 6
        # JacknifeAcrossCancerTypes('CCLE', DataType='Expression', MutType='KsKa')
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLEExpressionKsKa*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df[df['Coefficient'] == 'LogScore'] 
            df = df.merge(Complexes, right_on='Hugo', left_on='GeneName')
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'WithinCancer_TCGA' : # Rev Response #1 - Within cancer type GLMM box plot
        #   GetWithinCancerType(Dataset='TCGA', DataType='Protein')
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/WithinCancerType/TCGAExpressionKsKa*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df[df['Coefficient'] == 'LogScore'] 
            df = df.merge(Complexes, right_on='Hugo', left_on='GeneName')
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'WithinCancerGrouped_TCGA': # Rev Response #1 - Within cancer type GLMM heat map
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/WithinCancerType/TCGAExpressionKsKaByComplexGroup*'))
        NumSamples = GetTissueType('TCGA').groupby('type').size().reset_index()
        Metadata = AnnotateMutationalLoad(GetPointMutations('TCGA'),'KsKa')
        Metadata = Metadata.merge(GetTissueType('TCGA'), left_on='Barcode', right_on='Barcode').groupby('type')['MutLoad'].median().reset_index()
        Metadata = Metadata.merge(GetTissueType('TCGA').groupby('type').size().reset_index().rename(columns={0:'NumSamples'}))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['type'] = FileName.split('/')[11].replace('TCGAExpressionKsKaByComplexGroupOLSRegressionWithinCancerType','')
            df = df.merge(Metadata, left_on='type', right_on='type')
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'TMB_Shuffling_BetaDist': # Rev Response #2 - TMB shuffle + visalize p val dist
        shuff = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesTCGAShuffledTMB').assign(Group='Null')
        obs = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity').assign(Group='Observed')
        df = shuff.append(obs)
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['GeneName'] = df['GeneName'].str.split('_', expand=True)[0]
        return(ConvertPandasDFtoR(df))
        #df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['Pr...t..']), method = 'fdr')
        #GeneSet = df[(df['Estimate'] > 0) & (df['Adj.Pval'] < 0.05)]['GeneName']
        #gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(GeneSet)))
        #return(ConvertPandasDFtoR(gse[['term_name','source','p_value']]))
    elif FigureName == 'TMB_Shuffling_SigGeneCount':
        shuff = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesTCGAShuffledTMB').assign(Group='Null')
        obs = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity').assign(Group='Observed')
        df = shuff.append(obs)
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['GeneName'] = df['GeneName'].str.split('_', expand=True)[0]
        return(ConvertPandasDFtoR(df))


def GetFigure(Figure):
    SetUpPlottingPackages()
    if Figure == 'GlobalGSE_TCGA_Regression':  # Fig 1B
        ro.r.PlotCircularCORUMTCGA(GetFigureInput('GlobalGSE_TCGA_Regression'))
    elif Figure == 'GlobalGSE_TCGA_Regression_KEGG': # Fig 1C
        ro.r.PlotCircularKeggTCGA(GetFigureInput('GlobalGSE_TCGA_Regression'))
    elif Figure == 'AS_PSI': # Fig 2A
        # Raw data generated from fxn in GetExpLevelsForASEvents(ASType, PSI_Threshold) in Splicing.py
        ro.r.VisualizeAS(GetFigureInput('AS_PSI'))
    elif Figure == 'AS_Delta_PSI': # Fig 2B
        ro.r.PlotDeltaPSI(GetFigureInput('AS_Delta_PSI'))
    elif Figure == 'Groups_CCLEAndTCGA': # Fig 3
        ro.r.PlotRegCoefPerGroup(GetFigureInput('Groups_CCLEAndTCGA'))
    elif Figure == 'JacknifeshRNACCLE': # Fig 4A
        ro.r.PlotJacknifedshRNAAcrossGroups(GetFigureInput('JacknifeshRNACCLE'))
    elif Figure == 'JacknifeDrugsCCLE': # fig 4B
        ro.r.PlotJacknifedDrugsAcrossGroups(GetFigureInput('JacknifeDrugsCCLE'))
    elif Figure == 'AllDrugsByLoad': # Fig 5A
        ro.r.PlotFractionOfSigDrugsForAll(GetFigureInput('AllDrugsByLoad'))
    elif Figure == 'BootstrappedDrugs': # Fig 5B
        ro.r.PlotBootstrappedNegativelyAssociatedDrugs()
    elif Figure == 'Multicollinearity': # Sup Fig 1
        ro.r.PlotMultiCollinearityOfSNVsAndCNVs(GetFigureInput('Multicollinearity'))
    elif Figure == 'NumASEventsFiltered': # Sup Fig 2A
        ro.r.PlotNumASEventsFiltered(GetFigureInput('NumASEventsFiltered'))
    elif Figure == 'AS_PSI_NotFiltered':  # Sup Fig 2B
        ro.r.VisualizeAS(GetFigureInput('AS_PSI_NotFiltered'))
    elif Figure == 'AS_PSI': # Sup Fig 2C
        ro.r.VisualizeAS(GetFigureInput('AS_PSI'))
    elif Figure == 'PSI_Supplemental': # Sup Fig 3
        ro.r.VisualizeAllASThresholds() # 3A
        ro.r.VisualizeAllAS() # 3B
    elif Figure == 'JacknifeExpTCGA':  # Sup Fig 4
        ro.r.PlotJacknifedExpressionAcrossGroups(GetFigureInput('JacknifeExpTCGA'),'TCGA')
    elif Figure == 'CorrelationWithAge': # Sup Fig 5
        ro.r.PlotGLMMRegressionCoefficientsByAge(GetFigureInput('CorrelationWithAge'))
    elif Figure == 'JacknifeExpCCLE': # Sup Fig 6
        ro.r.PlotJacknifedExpressionAcrossGroups(GetFigureInput('JacknifeExpCCLE'), 'CCLE')
    elif Figure == 'RevResp_Within_TCGA_Exp': # Reviewer response - within cancer type
        ro.r.PlotWithinCancerTypeExpressionAcrossGroups(GetFigureInput('WithinCancer_TCGA'),'TCGA')
    elif Figure == 'RevResp_WithinCancerGrouped_TCGA': # Reviewer response - within cancer type
        ro.r.PlotCancerTypeInComplexesAcrossGroups(GetFigureInput('WithinCancerGrouped_TCGA'))
    elif Figure == 'RevResp_BetaDistAfterShuffle': # Reviewer response - p-val + beta dist in null and obs
        ro.r.VisualizeBetaAfterShuffle(GetFigureInput('TMB_Shuffling_BetaDist'))
    elif Figure == 'RevResp_AdjustedPValDist': # Reviewer response - adj and og pval dist
        ro.r.ComparePValueDistributions(GetFigureInput('TMB_Shuffling_SigGeneCount'))  
    elif Figure == 'RevResp_AdjustedPValGeneSetEnrichment': #Reviewer response - adj and g pval gse
        ro.r.PlotAdjustedPValGeneSetEnrich(GetFigureInput('TMB_Shuffling_SigGeneCount'))

        