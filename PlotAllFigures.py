
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
    rstats = importr('stats')
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
        pd.DataFrame({'Group': 'Immunoproteasome Core' , 'subgroup' : 'Proteasome', 'Hugo' : corum[corum['ComplexName'] == 'Immunoproteasome']['subunits(Gene name)'].str.split(';', expand=True).melt()['value'].str.upper() }),
        pd.DataFrame({'Group': '11S Regulatory Particle' , 'subgroup' : 'Proteasome', 'Hugo' :  ['PSME1','PSME2','PSME3'] }),
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
    InputDir = '/labs/ccurtis2/tilk/scripts/protein/Data/'
    if FigureNum == 'GlobalGSE_TCGA_Regression':
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ro.r.GetGeneSets('ExpressionTCGA')))
        gse = gse[gse['source'] == 'CORUM']
        return(ConvertPandasDFtoR(gse[['term_name','p_value']]))
    elif FigureNum == 'DGE_GSE_TCGA':
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ro.r.GetGeneSets('PanCancerDGE')))
        gse = gse[gse['source'] == 'CORUM']
        return(ConvertPandasDFtoR(gse[['term_name','p_value']]))
    elif FigureNum == 'GlobalGSE_CCLE_Regression':   
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment(ro.r.GetGeneSets('ExpressionCCLE')))
        gse = gse[gse['source'] == 'CORUM']
        return(ConvertPandasDFtoR(gse[['term_name','p_value']]))
    elif FigureNum == 'Groups_ByGene_TCGA':
        all = pd.read_csv(InputDir + 'Regression/ExpressionMixedEffectRegressionEstimatesKsKaTCGA')
        all['Coefficient'] = all['Unnamed: 0'].str.replace('\d+', '')
        all = all[all['Coefficient'] == 'LogScore']
        groups = GetGeneAnnotationsOfInterest()
        all = all.merge(groups, left_on='GeneName', right_on='Hugo')
        return(ConvertPandasDFtoR(all.astype(str)))   
    elif FigureNum == 'Groups_ByGene_CCLE':
        all = pd.read_csv(InputDir + 'Regression/ExpressionOLSRegressionEstimatesKsKaCCLE')
        all['Coefficient'] = all['Unnamed: 0'].str.replace('\d+', '')
        all = all[all['Coefficient'] == 'LogScore']
        groups = GetGeneAnnotationsOfInterest()
        all = all.merge(groups, left_on='GeneName', right_on='Hugo')
        return(ConvertPandasDFtoR(all.astype(str)))  
    elif FigureNum == 'AS':
        splicing = pd.read_csv(InputDir + 'Regression/ASMixedEffectRegressionEstimatesKsKaTCGA')
        splicing['Coefficient'] = splicing['Unnamed: 0'].str.replace('\d+', '')
        splicing = splicing[splicing['Coefficient'] == 'LogScore']
        splicing['GeneName'] = splicing['GeneName'].str.split('_', expand=True)[0]
        splicing['Adj.Pval'] = rstats.p_adjust(FloatVector(splicing['pVal']), method = 'fdr')
        GeneSet = splicing[(splicing['Estimate'] > 0) & (splicing['Adj.Pval'] < 0.05)]['GeneName'].unique().tolist()
        gse = ConvertRDataframetoPandas(ro.r.DoGeneSetEnrichment( ConvertPandasDFtoR(GeneSet)))
        gse[gse['source'] == 'REAC']
    elif FigureNum == 'Protein_CCLE':
        protein = pd.read_csv(InputDir + 'Regression/ProteinOLSRegressionEstimatesKsKaCCLE')
        protein['Coefficient'] = protein['Unnamed: 0'].str.replace('\d+', '')
        protein = protein[protein['Coefficient'] == 'LogScore']
        protein['Adj.Pval'] = rstats.p_adjust(FloatVector(protein['pVal']), method = 'fdr')
        protein = protein.merge(GetGeneAnnotationsOfInterest(), left_on='GeneName', right_on='Hugo')
        return(ConvertPandasDFtoR(protein.astype(str)))
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
    elif FigureNum == 'Chaperome_RNAi':
        all = pd.read_csv(InputDir + 'Regression/RNAiOLSRegressionEstimatesCCLE')
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/GeneSets/Human_Chaperome_TableS1A_PMID_29293508', sep='\t')
        all = all.merge(df[['Gene','Level2']], left_on= 'GeneName', right_on='Gene', how='left')
    elif FigureNum == 'All_RNAi':
        df = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/StandardizedRNAiOLSRegressionEstimatesKsKaCCLE')
        df['Coefficient'] = df['Unnamed: 0'].str.replace('\d+', '')
        df = df[df['Coefficient'] == 'NormalizedMutLoad']
        df[(df['adj.pval'] < 0.05) & (df['Estimate'] < 0 )]
   
   

def GetFigure(Figure):
    if Figure == 'Stability':
        foo=GetFigureInput('Stability')
        SetUpPlottingPackages()
        ro.r.PlotStabilityPerTumor(foo)
    elif Figure == 'Stability_DeltaG':
        foo=GetFigureInput('Stability_DeltaG')
        SetUpPlottingPackages()
        ro.r.PlotDeltaGPerTumor(foo)
    elif Figure == 'Groups_CCLE':
        ro.r.PlotRegCoefPerGroup(GetFigureInput('Groups_ByGene_CCLE'), 'CCLE')
    elif Figure == 'Protein_CCLE':
        ro.r.PlotRegCoefPerGroup(GetFigureInput('Protein_CCLE'), 'Protein_CCLE')
    elif Figure == 'Groups_TCGA':
        ro.r.PlotRegCoefPerGroup(GetFigureInput('Groups_ByGene_TCGA'), 'TCGA')
    elif Figure == 'GlobalGSE_TCGA_Regression':
        ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_TCGA_Regression'), 'TCGA_Regression')
    elif Figure == 'GlobalGSE_CCLE_Regression':
        ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_CCLE_Regression'), 'CCLE_Regression')
    elif Figure == 'DGE_GSE_TCGA':
        ro.r.PlotGlobalGSETCGA(GetFigureInput('DGE_GSE_TCGA'), 'TCGA_DGE')



SetUpPlottingPackages(); GetFigure(Figure='Groups_TCGA')

SetUpPlottingPackages(); GetFigure(Figure='Groups_CCLE')

SetUpPlottingPackages(); GetFigure(Figure='Protein_CCLE')




### GENE SET ENRICHMENT PLOTS -- FIG 1
SetUpPlottingPackages(); ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_TCGA_Regression'))
SetUpPlottingPackages(); ro.r.PlotGlobalGSETCGA(GetFigureInput('GlobalGSE_CCLE_Regression'))
SetUpPlottingPackages(); ro.r.PlotGlobalGSETCGA(GetFigureInput('DGE_GSE_TCGA'))








