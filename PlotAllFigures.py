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
from GetData import * 

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
        # Aggregstes all OLS regression with purity as explanatory
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/WithinCancerType/TCGAExpressionKsKaMixedEffectRegressionWithinCancer*'))
        NumSamples = GetTissueType('TCGA').groupby('type').size().reset_index()
        Metadata = AnnotateMutationalLoad(GetPointMutations('TCGA'),'KsKa')
        Metadata = Metadata.merge(GetTissueType('TCGA'), left_on='Barcode', right_on='Barcode').groupby('type')['MutLoad'].median().reset_index()
        Metadata = Metadata.merge(GetTissueType('TCGA').groupby('type').size().reset_index().rename(columns={0:'NumSamples'}))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df['type'] = FileName.split('/')[11].replace('TCGAExpressionKsKaMixedEffectRegressionWithinCancerType','')
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df[df['Coefficient'] == 'LogScore'] 
            df = df.merge(Metadata, left_on='type', right_on='type')
            Out = Out.append(df)
        Out = Out.merge(Complexes, right_on='Hugo', left_on='GeneName')
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
    elif FigureName == 'TMB_Shuffling_SigGeneCount':# Rev Response #2 - TMB shuffle + gse
        shuff = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesTCGAShuffledTMB').assign(Group='Null')
        obs = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity').assign(Group='Observed')
        df = shuff.append(obs)
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['GeneName'] = df['GeneName'].str.split('_', expand=True)[0]
        return(ConvertPandasDFtoR(df))
    elif FigureName == 'ResidualFit_Expression_TCGA': # Rev response - WIP raw residual visualization
        model = pd.read_csv(InputDir + 'Regression/Diagnostics/CCLEExpressionKsKaModelDiagnosticResidualVsFittedVals')
        df = model.merge(Complexes, right_on='Hugo', left_on='GeneName', how='left' )
    elif FigureName == 'Protein_CCLE_Subset': # Rev response - protein expression for subset in CCLE (within cancer type + all)
        Out=pd.DataFrame()
        tissue = GetTissueType('CCLE').groupby(['type']).size().reset_index().rename(columns={0:'NumSamples'})
        df = pd.read_csv(InputDir + 'Regression/ProteinOLSRegressionEstimatesKsKaCCLE')        
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        all = df.merge(Complexes, right_on='Hugo', left_on='GeneName').assign(type='All').assign(NumSamples=tissue['NumSamples'].sum())
        all = all[['adj.r.squared', 'sigma', 'Estimate', 'Std..Error', 't.value', 'Pr...t..', 'pVal', 'GeneName', 'type', 'Group', 'subgroup', 'Hugo','NumSamples', 'Coefficient']]
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/WithinCancerType/CCLEProteinKsKaOLSRegressionWithinCa*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName).rename(columns={'CancerType':'type'})
            df = df.merge(Complexes, left_on='GeneName', right_on='Hugo')
            df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
            df = df.merge(tissue,left_on= 'type', right_on = 'type', how='left')
            df = df[df['Coefficient'] == 'LogScore']
            df = df[['adj.r.squared', 'sigma', 'Estimate', 'Std..Error', 't.value', 'Pr...t..', 'pVal', 'GeneName', 'type', 'Group', 'subgroup', 'Hugo', 'NumSamples','Coefficient']]
            Out = Out.append(df)
        Out = all.append(Out)
        return(ConvertPandasDFtoR(Out))
    elif FigureName == 'Protein_Vs_Exp_CCLE_IndividualGene':
        pro = pd.read_csv(InputDir + 'Regression/MedianProteinOLSRegressionEstimatesKsKaCCLE')
        pro['Coefficient'] = pro['Unnamed: 0'].str.replace('\d+', '')
        pro = pro[pro['Coefficient'] == 'LogScore']
        pro['AdjPval'] = rstats.p_adjust(FloatVector(pro['pVal']), method = 'fdr')   
        pro['DataType'] = 'Protein'
        exp = pd.read_csv(InputDir + 'Regression/ExpressionOLSRegressionEstimatesKsKaCCLE')  
        exp['Coefficient'] = exp['Unnamed: 0'].str.replace('\d+', '')
        exp = exp[exp['Coefficient'] == 'LogScore']
        exp['AdjPval'] = rstats.p_adjust(FloatVector(exp['pVal']), method = 'fdr')   
        exp['DataType'] = 'Expression'
        all = pro.append(exp)
        quantile = all.groupby(['DataType'])['Estimate'].quantile([0,0.1,0.5,0.9,1]).reset_index().rename(columns={'level_1':'Group'}).assign(
            subgroup='Quantile').assign(AdjPval = 0).assign(GeneName='')
        all = all.merge(Complexes, left_on='GeneName', right_on='Hugo')
        all = pd.concat([all[['DataType','Group','Estimate','subgroup','AdjPval','GeneName']], quantile])
        return(ConvertPandasDFtoR(all.astype(str)))
    elif FigureName == 'Protein_Vs_Exp_TCGA_BRCA_IndividualGene':
        #pro = pd.read_csv(InputDir + '/Regression/ProteinOLSRegressionWithPurityBRCAEstimatesKsKaTCGA') # TCGA BRCA protein only
        #pro = pd.read_csv(InputDir + '/Regression/ProteinOLSRegressionWithPurityBRCAandOVEstimatesKsKaTCGA') # TCGA BRCA and OV protein only
        pro = pd.read_csv(InputDir + '/Regression/ProteinMixedEffectRegressionWithoutPurityEstimatesKsKaCPTAC') # All avaiable CPTAC from cbioportal
        pro['Coefficient'] = pro['Unnamed: 0'].str.replace('\d+', '')
        pro = pro[pro['Coefficient'] == 'LogScore']
        pro['AdjPval'] = rstats.p_adjust(FloatVector(pro['Pr...t..']), method = 'fdr')   
        pro['DataType'] = 'Protein'
        exp = pd.read_csv(InputDir + 'Regression/ERChapAdded/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')  
        exp['Coefficient'] = exp['Unnamed: 0'].str.replace('\d+', '')
        exp = exp[exp['Coefficient'] == 'LogScore']
        exp['AdjPval'] = rstats.p_adjust(FloatVector(exp['Pr...t..']), method = 'fdr')   
        exp['DataType'] = 'Expression'
        all = pro.append(exp)
        quantile = all.groupby(['DataType'])['Estimate'].quantile([0,0.1,0.5,0.9,1]).reset_index().rename(columns={'level_1':'Group'}).assign(
            subgroup='Quantile').assign(AdjPval = 0).assign(GeneName='')
        all = all.merge(Complexes, left_on='GeneName', right_on='Hugo')
        all = pd.concat([all[['DataType','Group','Estimate','subgroup','AdjPval','GeneName']], quantile])
        return(ConvertPandasDFtoR(all.astype(str)))
    elif FigureName == 'AdjPValueAllCorum': # adjusted p-value between all complexes in CORUM in shuffled vs observed - RR
        # shuff = pd.DataFrame() # empty df to append all shuffled permutations to 
        # ShuffFiles = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Shuffling/*'))
        # Rep = 0 #
        # for FileName in ShuffFiles: # merge all permutations/shuffled TMB estimates into one df 
        #     Rep = Rep + 1
        #     shuff = shuff.append(pd.read_csv(FileName).assign(Group = 'Null_' + str(Rep)))
        # obs = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity').assign(Group='Observed')
        # df = shuff.append(obs) #  add shuffled to the observed non shuffled regression estimate
        # df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        # df = df[df['Coefficient'] == 'LogScore']
        # Out = pd.DataFrame() # empty df to append final results
        # for gene in df['GeneName'].unique():
        #     SingleGene = df[df['GeneName'] == gene]
        #     NullSingleGene = SingleGene[SingleGene['Group'] != 'Observed']
        #     ObsEstimate = SingleGene[SingleGene['Group'] == 'Observed']['Estimate'].tolist()[0] # beta coefficient for gene in observed/non-shuffled data
        #     if ObsEstimate < 0:
        #         AdjPval = len(NullSingleGene[NullSingleGene['Estimate'] < ObsEstimate])/len(NullSingleGene)
        #     else: # observed estimate for gene is positive 
        #         AdjPval = len(NullSingleGene[NullSingleGene['Estimate'] > ObsEstimate])/len(NullSingleGene)
        #     Out = Out.append( SingleGene[SingleGene['Group'] == 'Observed'].assign(AdjPval = AdjPval).reset_index())
        # corum = pd.read_csv(os.getcwd() + '/GeneSets/coreComplexes.txt.zip', sep='\t')
        # corum = corum.set_index(['ComplexName'])['subunits(Gene name)'].str.split(';', expand=True).fillna(value=np.nan).reset_index()
        # corum = pd.melt(corum, id_vars=['ComplexName']).dropna().rename(columns={'value':'GeneName'})
        # corum['GeneName'] = corum['GeneName'].str.upper()
        # Out = Out.merge(corum[['GeneName','ComplexName']], right_on='GeneName', left_on='GeneName', how='left')#.rename(columns={'Estimate':'Null_Estimate'})
        #Out.to_csv('/home/tilk/ShuffledByGene')
        Out = pd.read_csv('/home/tilk/ShuffledByGene')
        return(ConvertPandasDFtoR(Out.dropna()))
    elif FigureName == 'CCLE_Correlation_RNA_and_Protein':
        pro = pd.read_csv(InputDir + 'Regression/MedianProteinOLSRegressionEstimatesKsKaCCLE').rename(columns={'Estimate':'Protein_Estimate'})
        pro['Coefficient'] = pro['Unnamed: 0'].str.replace('\d+', '')
        pro = pro[pro['Coefficient'] == 'LogScore']
        pro['AdjPval'] = rstats.p_adjust(FloatVector(pro['pVal']), method = 'fdr')   
        pro['DataType'] = 'Protein'
        exp = pd.read_csv(InputDir + 'Regression/ExpressionOLSRegressionEstimatesKsKaCCLE').rename(columns={'Estimate':'Exp_Estimate'})
        exp['Coefficient'] = exp['Unnamed: 0'].str.replace('\d+', '')
        exp = exp[exp['Coefficient'] == 'LogScore']
        all = pro[['Protein_Estimate','GeneName']].merge(exp[['Exp_Estimate','GeneName']], left_on='GeneName', right_on='GeneName')
        all = all.merge(Complexes, left_on='GeneName', right_on='Hugo')
        return(ConvertPandasDFtoR(all.astype(str)))
    elif FigureName == 'RawResidualsInGenesOfInterest_TCGA':
        df = pd.read_csv(InputDir + 'Regression/Diagnostics/TCGAExpressionKsKaModelDiagnosticResidualVsFittedVals')
    elif FigureName == 'VarianceExpressionByCancerType':
        from scipy.stats import zscore
        df = GetExpressionData('TCGA')
        df = df.merge(GetTissueType('TCGA'), left_on='Barcode', right_on='Barcode')
        df['NormalizedValue'] = df.groupby(['GeneName']).Value.transform(lambda x : zscore(x,ddof=1))
        df = df.groupby(['type','GeneName'])['Value'].var().reset_index()
        df = df.merge(Complexes, left_on='GeneName', right_on='Hugo', how='left')
        df.to_csv(InputDir + 'Regression/TCGAVarianceInGeneExpression')



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
    elif Figure == 'Protein_CCLE_Subset': # Rev response - protein in complexes in CCLE
        ro.r.PlotProteinAcrossCancerTypesInGroups(GetFigureInput('Protein_CCLE_Subset'))
    elif Figure == 'Protein_Vs_Exp_CCLE_IndividualGene': # Rev response - individual gene protein vs exp
        ro.r.ProteinVsExpressionByIndividualGene(GetFigureInput('Protein_Vs_Exp_CCLE_IndividualGene'),'CCLE')
    elif Figure == 'Protein_Vs_Exp_TCGA_BRCA_IndividualGene': # Rev response - individual gene protein vs exp
        ro.r.ProteinVsExpressionByIndividualGene(GetFigureInput('Protein_Vs_Exp_TCGA_BRCA_IndividualGene'),'TCGA')
    elif Figure == 'AdjPValueAllCorum': # Rev response -- p-value for all corum, no gse
        ro.r.PlotAllCorumNullvsObs(GetFigureInput('AdjPValueAllCorum'))
    elif Figure == 'CCLE_Correlation_RNA_and_Protein': # RR - correlation coefficients between rna and protein
        ro.r.PlotCorrelationCoefficientsBetweenExpAndProt(GetFigureInput('CCLE_Correlation_RNA_and_Protein'))


#GetFigureInput('AdjPValueAllCorum')

GetFigureInput('VarianceExpressionByCancerType')


# ro.r.PlotCancerTypeInComplexesAcrossGroups(GetFigureInput('WithinCancerGrouped_TCGA'))