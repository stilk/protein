
import pandas as pd
import os
import numpy as np
from pandas.core.frame import DataFrame
import rpy2.robjects as ro 
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.vectors import FloatVector
from GetData import * 
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
        pd.DataFrame({'Group': 'Cytoplasmic Ribosomes','subgroup' : 'Translation', 'Hugo' : corum[corum['ComplexName'].str.contains('subunit, cytoplasmic')]['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].drop_duplicates().tolist()}),
        pd.DataFrame({'Group': 'Mitochondrial Ribosomes', 'subgroup' : 'Translation', 'Hugo' : corum[corum['ComplexName'].str.contains('subunit, mitochondrial')]['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].drop_duplicates().tolist()}),
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
    if FigureNum == 'GlobalGSE_TCGA_Regression':
        df = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity')
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['GeneName'] = df['GeneName'].str.split('_', expand=True)[0]
        df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['Pr...t..']), method = 'fdr')
        GeneSet = df[(df['Estimate'] > 0) & (df['Adj.Pval'] < 0.05)]['GeneName']
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(GeneSet)))
        #gse = gse[gse['source'] == 'CORUM']
        return(ConvertPandasDFtoR(gse[['term_name','source','p_value']]))
    elif FigureNum == 'DGE_GSE_TCGA':
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ro.r.GetGeneSets('PanCancerMatchedDGE')))
        gse = gse[gse['source'] == 'CORUM']
        return(ConvertPandasDFtoR(gse[['term_name','p_value']]))
    elif FigureNum == 'GlobalGSE_CCLE_Regression':   
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ro.r.GetGeneSets('ExpressionCCLE')))
        gse = gse[gse['source'] == 'CORUM']
        return(ConvertPandasDFtoR(gse[['term_name','p_value']]))
    elif FigureNum == 'Groups_CCLEAndTCGA':
        tcga = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGA').assign(Dataset='TCGA')
        ccle = pd.read_csv(InputDir + 'Regression/ExpressionOLSRegressionEstimatesKsKaCCLE').assign(Dataset='CCLE')
        all = pd.concat([tcga[['Pr...t..','adj.pval','Estimate','GeneName','Unnamed: 0','Dataset']], 
                    ccle[['Pr...t..','adj.pval','Estimate','GeneName','Unnamed: 0','Dataset']]])
        all['Coefficient'] = all['Unnamed: 0'].str.replace('\d+', '')
        all = all[all['Coefficient'] == 'LogScore']
        groups = GetGeneAnnotationsOfInterest()
        all = all.merge(groups, left_on='GeneName', right_on='Hugo')
        return(ConvertPandasDFtoR(all.astype(str))) 
    elif FigureNum == 'GlobalTranscription_TCGA': 
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_AvgGeneExpForAllCancerTypes')
        Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')[['Barcode','MutLoad']]
        df = df.merge(Samples, left_on='Barcode', right_on='Barcode', how='left')
        df['Bin'] = pd.cut(df['MutLoad'], bins=[0,100,1000,50000], labels=['0 - 100','100 - 1,000','1,0000 - >10,000'])
        df['Housekeeping'] =df['GeneName'].isin(pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/housekeeping_genes.txt.gz', sep='   ')['HUGO'])
        df['Essential'] = df['GeneName'].isin(pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/essential_genes.txt.gz', header=None)[0])
        df['Drivers'] = df['GeneName'].isin(pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/drivers_Bailey2018.txt', header=None)[0])
        over = pd.concat([df.query('OverExpressedTranscript == True').groupby(['Bin','Housekeeping','Essential','Drivers']).size(), 
                    df.query('OverExpressedTranscript == False').groupby(['Bin','Housekeeping','Essential','Drivers']).size()], 
            axis=1).reset_index().rename(columns={0:'OverExpressed',1:'NotOverExpressed'})
        over['FractionOverExpressedGenes'] = over['OverExpressed']/(over['OverExpressed'] + over['NotOverExpressed'])
        under = pd.concat([df.query('DepletedTranscript == True').groupby(['Bin','Housekeeping','Essential','Drivers']).size(), 
                df.query('DepletedTranscript == False').groupby(['Bin','Housekeeping','Essential','Drivers']).size()], 
            axis=1).reset_index().rename(columns={0:'UnderExpressed',1:'NotUnderExpressed'})
        under['FractionUnderExpressedGenes'] = under['UnderExpressed']/(under['UnderExpressed'] + under['NotUnderExpressed'])
        out = under.merge(over, left_on=['Bin','Housekeeping','Essential','Drivers'],
            right_on=['Bin','Housekeeping','Essential','Drivers'])
        out.to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_OverAndExpressionCountsFinerBinsAndGeneSets')
    elif FigureNum == 'GlobalTranscription_CCLE': 
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/CCLE_AvgGeneExpForAllCancerTypes')
        Samples = AnnotateMutationalLoad(GetMutations('CCLE'), 'KsKa')
        df = df.merge(Samples, left_on='Barcode', right_on='Barcode', how='left')
        df['Bin'] = pd.cut(df['MutLoad'], bins=[0,100,1000,50000], labels=['0 - 100','100 - 1,000','1,0000 - >10,000'])
        df['Housekeeping'] =df['GeneName'].isin(pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/housekeeping_genes.txt.gz', sep='   ')['HUGO'])
        df['Essential'] = df['GeneName'].isin(pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/essential_genes.txt.gz', header=None)[0])
        df['Drivers'] = df['GeneName'].isin(pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/drivers_Bailey2018.txt', header=None)[0])
        over = pd.concat([df.query('OverExpressedTranscript == True').groupby(['Bin','Housekeeping','Essential','Drivers']).size(), 
                    df.query('OverExpressedTranscript == False').groupby(['Bin','Housekeeping','Essential','Drivers']).size()], 
            axis=1).reset_index().rename(columns={0:'OverExpressed',1:'NotOverExpressed'})
        over['FractionOverExpressedGenes'] = over['OverExpressed']/(over['OverExpressed'] + over['NotOverExpressed'])
        under = pd.concat([df.query('DepletedTranscript == True').groupby(['Bin','Housekeeping','Essential','Drivers']).size(), 
                df.query('DepletedTranscript == False').groupby(['Bin','Housekeeping','Essential','Drivers']).size()], 
            axis=1).reset_index().rename(columns={0:'UnderExpressed',1:'NotUnderExpressed'})
        under['FractionUnderExpressedGenes'] = under['UnderExpressed']/(under['UnderExpressed'] + under['NotUnderExpressed'])
        out = under.merge(over, left_on=['Bin','Housekeeping','Essential','Drivers'],
            right_on=['Bin','Housekeeping','Essential','Drivers'])
        out.to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/CCLE_OverAndExpressionCountsFinerBinsAndGeneSets') 
    elif FigureNum == 'AS_Regression':
        splicing = pd.read_csv(InputDir + 'Regression/ASMixedEffectRegressionEstimatesKsKaTCGA')
        splicing['Coefficient'] = splicing['Unnamed: 0'].str.replace('\d+', '')
        splicing = splicing[splicing['Coefficient'] == 'LogScore']
        splicing['GeneName'] = splicing['GeneName'].str.split('_', expand=True)[0]
        splicing['Adj.Pval'] = rstats.p_adjust(FloatVector(splicing['pVal']), method = 'fdr')
        GeneSet = splicing[(splicing['Estimate'] > 0) & (splicing['Adj.Pval'] < 0.05)]['GeneName'].unique().tolist()
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(GeneSet)))
        gse[gse['source'] == 'REAC']
    elif FigureNum == 'AS_Increase_Load':
        OutDir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
        ASType='ES'
        df = GetAlternativeSplicingData(AS=ASType).stack().reset_index().rename(columns={'level_6':'Barcode',0:'PSI'})
        df = df[~df['Barcode'].str.contains('Norm')]
        drivers = pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/drivers_Bailey2018.txt', header=None)[0]
        df['driver'] = df['symbol'].isin(drivers)
        Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
        Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
        df = df.merge(Samples, left_on='Barcode', right_on='Barcode', how='left')
        df.groupby(['Barcode','MutLoad','driver'])['PSI'].mean().reset_index().to_csv(OutDir + 'TCGA_' + ASType + '_Avg_PSI_PerTumor_ByDrivers')
        df.groupby(['Barcode','MutLoad'])['PSI'].mean().reset_index().to_csv(OutDir + 'TCGA_'+ ASType +'_Avg_PSI_PerTumor')
        df = df.dropna()
        df = pd.concat( [ df[df['PSI'] <= 0.1].assign(Group='Not_Retained'), df[df['PSI'] >= 0.9].assign(Group='Retained')])
        df.groupby(['Barcode','MutLoad','Group']).size().unstack().reset_index().to_csv(OutDir + 'TCGA_'+ ASType +'_Counts_ThresholdByPSI_PerTumor')
        df.groupby(['Barcode','MutLoad','driver','Group']).size().unstack().reset_index().to_csv(OutDir + 'TCGA_'+ ASType +'_Counts_ThresholdByPSI_PerTumorAndDrivers')
    elif FigureNum == 'AS_Increase_Load_TSLevel':
        foo = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_AS_RI_AvgExpEffect', low_memory=False)
        df = pd.concat( [ foo[foo['PSI'] <= 0.1].assign(Group='Not_Retained'), foo[foo['PSI'] >= 0.9].assign(Group='Retained')])
        drivers = pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/drivers_Bailey2018.txt', header=None)[0]
        df['driver'] = df['GeneName'].isin(drivers)
        Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
        Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
        df = df.merge(Samples, left_on='Barcode', right_on='Barcode', how='left')
        OutDir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'; ASType='RI'
        df.groupby(['Barcode','MutLoad','Group','DepletedTranscript']).size().unstack().reset_index().to_csv(OutDir + 'TCGA_'+ ASType +'_Counts_ThresholdByPSI_PerTumor')
        df.groupby(['Barcode','MutLoad','driver','Group','DepletedTranscript']).size().unstack().reset_index().to_csv(OutDir + 'TCGA_'+ ASType +'_Counts_ThresholdByPSI_PerTumorAndDrivers')
    #***************
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
    elif FigureNum == 'AllByGene_TCGA':
        all = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGA')
        all['Coefficient'] = all['Unnamed: 0'].str.replace('\d+', '')
        all = all[all['Coefficient'] == 'LogScore']
        all = all[(all['adj.pval'] < 0.05) & (all['Estimate'] > 0)]
        all['pval_rank'] =all['adj.pval'].rank(method='min')
        corum = pd.read_csv('/labs/ccurtis2/tilk/02_DNDS/separateDatabaseFiles/CORUM/coreComplexes.txt.zip', sep='\t')
        PosGSE = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(all['GeneName'])))
        GenesInCorumEnrichment = corum[corum['ComplexName'].isin(PosGSE['term_name'])].set_index(['ComplexName'])['subunits(Gene name)'].str.split(
            ';', expand=True).reset_index().melt(id_vars='ComplexName').rename(columns={'value':"GeneName"})   
        out = all.merge(GenesInCorumEnrichment, left_on='GeneName', right_on='GeneName', how='left').replace(np.nan,0)
        out['NegLog10PVal'] = -np.log10(out['adj.pval'])
        out = out[['Estimate','ComplexName','NegLog10PVal','pval_rank','GeneName']]
        return(ConvertPandasDFtoR(out.astype(str)))
    elif FigureNum == 'Drug_CCLE':
        df = GetPerDrugRegression()
        df['Coefficient'] = df['Coefficient'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        return(ConvertPandasDFtoR(df))
    elif FigureNum == 'Protein_CCLE': 
        protein = pd.read_csv(InputDir + 'Regression/ProteinOLSRegressionEstimatesKsKaCCLE')
        protein['Coefficient'] = protein['Unnamed: 0'].str.replace('\d+', '')
        protein = protein[protein['Coefficient'] == 'LogScore']
        protein['Adj.Pval'] = rstats.p_adjust(FloatVector(protein['Pr...t..']), method = 'fdr')
        protein = protein.merge(GetGeneAnnotationsOfInterest(), left_on='GeneName', right_on='Hugo')
        return(ConvertPandasDFtoR(protein.astype(str)))
    elif FigureNum == 'Protein_Exp':
        protein = GetProteinExpressionData('CCLE')
        protein = protein.merge(GetGeneAnnotationsOfInterest(), left_on='GeneName', right_on='Hugo')
        Muts = AnnotateMutationalLoad(GetMutations('CCLE'), MutType='KsKa')
        Groups = pd.concat([Muts[Muts['MutLoad'] >= 2000].assign(MutGroup = 'High'), Muts[Muts['MutLoad'] <= 200].assign(MutGroup = 'Low')])
        out = protein.merge(Groups, left_on = 'Barcode', right_on = 'Barcode')
        return(ConvertPandasDFtoR(out.astype(str)))
    elif FigureNum == '1B': # TCGA Regression for nonsyn/polyphen/HE
        df = pd.concat([pd.read_csv(InputDir + '/Regression/TCGA_Expression_Chaperome_Nonsynonymous.txt').assign(Type='Nonsynonymous'),
                        pd.read_csv(InputDir + '/Regression/TCGA_Expression_Chaperome_Polyphen.txt').assign(Type='Damaging Nonsynonymous'),
                        pd.read_csv(InputDir + '/Regression/TCGA_Expression_Chaperome_SNVs.txt').assign(Type='All SNVs')])
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
    elif FigureNum == '1B_v2': # Grouped TCGA Regression for nonsyn/polyphen/HE
        df = pd.concat([pd.read_csv(InputDir + '/Regression/GroupedRegression_TCGA_Expression_Chaperome_Nonsynonymous.txt').assign(Type='Nonsynonymous'),
                        pd.read_csv(InputDir + '/Regression/GroupedRegression_TCGA_Expression_Chaperome_Polyphen.txt').assign(Type='Damaging Nonsynonymous'),
                        pd.read_csv(InputDir + '/Regression/GroupedRegression_TCGA_Expression_Chaperome_SNVs.txt').assign(Type='All SNVs')])
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
    elif FigureNum == '1B_QQ':
        df = pd.read_csv(InputDir + '/Regression/ExpressionMixedEffectRegressionEstimatesTCGA')
    elif FigureNum == 'Stability':
        mut =  AnnotateMutationalLoad(GetMutations(Dataset='TCGA'), MutType='SNV')
        stability = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/MutFunc/TCGA_ProteinStability_HomologyModel' )
        drivers = pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/drivers_Bailey2018.txt', header=None)[0]
        stability['Driver'] = stability['Hugo_Symbol'].isin(drivers)
        stability = stability[['MutID','ENST','Tumor_Sample_Barcode','Hugo_Symbol', 'Gene','Variant_Classification','ddG','Driver']].dropna()
        stability['Barcode'] = stability['Tumor_Sample_Barcode'].str[0:15]
        stability = stability.merge(mut, left_on='Barcode', right_on='Barcode', how='left')
        stability['DestabilizingMut'] = ((stability['ddG'] >= 1) * 1).map({1:'Unstable', 0:'Stable'})
        stability = stability.groupby(['Tumor_Sample_Barcode','MutLoad','Driver'])['DestabilizingMut'].value_counts().unstack().replace(np.nan,0).reset_index()
        return(ConvertPandasDFtoR(stability))
    elif FigureNum == 'Stability_DeltaG':
        mut =  AnnotateMutationalLoad(GetMutations(Dataset='TCGA'), MutType='SNV')
        stability = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/MutFunc/TCGA_ProteinStability_HomologyModel' )
        drivers = pd.read_csv('/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/drivers_Bailey2018.txt', header=None)[0]
        stability['Driver'] = stability['Hugo_Symbol'].isin(drivers)
        stability = stability[['MutID','ENST','Tumor_Sample_Barcode','Hugo_Symbol', 'Gene','Variant_Classification','ddG','Driver']].dropna()
        stability['Barcode'] = stability['Tumor_Sample_Barcode'].str[0:15]
        stability = stability.merge(mut, left_on='Barcode', right_on='Barcode', how='left')
        stability = stability.groupby(['Tumor_Sample_Barcode','MutLoad','Driver'])['ddG'].mean().reset_index()
        return(ConvertPandasDFtoR(stability))
        #return(stability)
        #stability['MutBin'] = pd.cut(stability['MutLoad'], bins=[0,10,1000,50000])
        # foo = stability.groupby(['MutBin'])[['Stable','Unstable']].sum().reset_index()
        # foo['Unstable']/(foo['Stable'] + foo['Unstable'])
        # stability['Unstable_Fraction'] = stability['Unstable']/(stability['Unstable'] + stability['Stable'])
    elif FigureNum == 'Grouped_RNAi':
        df = GetshRNARegression()
        df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['Pr...t..']), method = 'fdr')
        return(ConvertPandasDFtoR(df))
    elif FigureNum == 'All_RNAi':
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/RNAiOLSRegressionEstimatesNoNormKsKaCCLE')
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'LogScore']
        df['Adj.Pval'] = rstats.p_adjust(FloatVector(df['pVal']), method = 'fdr')
        test= df[(df['Adj.Pval'] < 0.05) & (df['Estimate'] < 0 )]['GeneName']
        test2=ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ConvertPandasDFtoR(test)))
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
    elif FigureNum == 'JacknifeDrugsCCLE':
        DrugsOfInterest = GetDrugResponseData(Screen='primary', AllDrugs=False)[['name','Group','subgroup']].drop_duplicates()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLEDrug*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df = df[df['Coefficient'] == 'LogScore'] 
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureNum == 'JacknifeshRNACCLE':
        Complexes = GetSubsetOfGeneAnnotationsOfInterest()
        Out = pd.DataFrame()
        ListOfCancerTypes = glob.glob(os.path.join(os.getcwd() + '/Data/Regression/Jacknife/CCLERNAi*'))
        for FileName in ListOfCancerTypes:
            df = pd.read_csv(FileName)
            df = df[df['Coefficient'] == 'LogScore'] 
            Out = Out.append(df)
        return(ConvertPandasDFtoR(Out))
    elif FigureNum == 'AllDrugsByLoad':
        df = pd.read_csv(os.getcwd() + '/Data/Regression/AllDrugsOLSRegressionEstimatesKsKaCCLE')
        df = df[df['Coefficient'] == 'LogScore'] 
        # Freq = df['subgroup'].value_counts().reset_index()
        # Freq['PlottingGroup'] = np.where(Freq['subgroup'] < 4,'Other', Freq['index'])
        # df = df.merge(Freq[['index','PlottingGroup']], left_on='subgroup',right_on='index')
        # df[df['pVal'] < 0.05]['PlottingGroup'].unique()
        # # CountsOfSubgroups = df['subgroup'].value_counts().reset_index()
        # CountsOfSubgroups = CountsOfSubgroups[CountsOfSubgroups['subgroup'] > 5]
        # df = df[df['subgroup'].isin(CountsOfSubgroups['index'])]
        return(ConvertPandasDFtoR(df))







def GetFigure(Figure):
    if Figure == 'Stability':
        foo=GetFigureInput('Stability')
        SetUpPlottingPackages()
        ro.r.PlotStabilityPerTumor(foo)
    elif Figure == 'Stability_DeltaG':
        foo=GetFigureInput('Stability_DeltaG')
        SetUpPlottingPackages()
        ro.r.PlotDeltaGPerTumor(foo)
    elif Figure == 'Groups_CCLEAndTCGA':
        all = GetFigureInput('Groups_CCLEAndTCGA')
        SetUpPlottingPackages(); ro.r.PlotRegCoefPerGroup(all)
    elif Figure == 'Protein_CCLE':
        ro.r.PlotRegCoefPerGroup(GetFigureInput('Protein_CCLE'), 'Protein_CCLE')
    elif Figure == 'Protein_Exp':
        ro.r.PlotProtein(GetFigureInput('Protein_Exp'))
    elif FigureNum == 'AS_Delta_PSI':
        out = GetFigureInput('AS_Delta_PSI')
        SetUpPlottingPackages(); ro.r.PlotDeltaPSI(out)
    elif FIgure == 'RNAi_CCLE':
        foo = GetFigureInput('Grouped_RNAi')
        SetUpPlottingPackages(); ro.r.PlotshRNA(foo)
    elif FigureNum == 'Drug_CCLE':
        ro.r.PlotDrug(GetFigureInput('Drug_CCLE'))
    elif Figure == 'AllIndividualGene_TCGA':
        ro.r.PlotRefCoefAllGenes(GetFigureInput('AllByGene_TCGA'))
        foo=GetFigureInput('AllByGene_TCGA')
        SetUpPlottingPackages(); ro.r.PlotRefCoefAllGenes(foo)
    elif Figure == 'GlobalGSE_TCGA_Regression':
        ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_TCGA_Regression'), 'TCGA_Regression')
    elif Figure == 'GlobalGSE_TCGA_Regression_KEGG':
        foo = GetFigureInput('GlobalGSE_TCGA_Regression')
        SetUpPlottingPackages(); ro.r.PlotGlobalKeggTCGA(foo)
    elif Figure == 'GlobalGSE_CCLE_Regression':
        ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_CCLE_Regression'), 'CCLE_Regression')
    elif Figure == 'DGE_GSE_TCGA':
        ro.r.PlotGlobalGSETCGA(GetFigureInput('DGE_GSE_TCGA'), 'TCGA_DGE')
    elif Figure == 'OverAndUnderGlobal':
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_OverAndExpressionCountsFinerBinsAndGeneSets')
        df['Drivers'] = df['Drivers'].map({False:'Passengers',True:'Drivers'})
        df['Essential'] = df['Essential'].map({False:'Not Essential',True:'Essential'})
        df['Housekeeping'] = df['Housekeeping'].map({False:'Not Housekeeping',True:'Housekeeping'})
        df = pd.melt(df, id_vars=['Bin','UnderExpressed', 'NotUnderExpressed', 'OverExpressed', 'NotOverExpressed'], 
            value_vars=['Housekeeping','Essential','Drivers'])
        df = df.groupby(['value','Bin'])[['UnderExpressed', 'NotUnderExpressed', 'OverExpressed','NotOverExpressed']].sum().reset_index()
        SetUpPlottingPackages(); ro.r.PlotGlobalDownAndUpregulation(ConvertPandasDFtoR(df))
    elif Figure == 'JacknifeExpCCLE':
        df = GetFigureInput('JacknifeExpCCLE')
        SetUpPlottingPackages(); ro.r.PlotJacknifedExpressionAcrossGroups(df, 'CCLE')
    elif Figure == 'JacknifeExpTCGA':
        df = GetFigureInput('JacknifeExpTCGA')
        SetUpPlottingPackages(); ro.r.PlotJacknifedExpressionAcrossGroups(df,'TCGA')
    elif Figure == 'JacknifeDrugsCCLE':
        df = GetFigureInput('JacknifeDrugsCCLE')
        SetUpPlottingPackages(); ro.r.PlotJacknifedDrugsAcrossGroups(df)
    elif Figure == 'JacknifeshRNACCLE':
        df = GetFigureInput('JacknifeshRNACCLE')
        SetUpPlottingPackages(); ro.r.PlotJacknifedshRNAAcrossGroups(df)
    elif Figure == 'AllDrugsByLoad':
        df = GetFigureInput('AllDrugsByLoad')
        SetUpPlottingPackages(); ro.r.PlotRegressionForAllDrugs(df)




        
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



