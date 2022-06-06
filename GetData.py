'''
This script reads in different data types (e.g., mutation calls, expression values, etc.) and metadata (e.g., tissue type) for both TCGA and CCLE.
'''

import pandas as pd
import numpy as np
import glob
import os
from scipy import stats
import scipy
import subprocess



def GetExpressionData(Dataset, GeneSet=[], DataDir = os.getcwd() + "/Data/Raw/Expression/"):
    '''
    Returns expression values for different datases of a dataframe with 3 columns: Barcode, GeneName, Value.
	Args: 	@Dataset = a string of which dataset to read values from ('CCLE', or 'TCGA')
    		@GeneSet = optional argument of a list of gene names to only read expression calls for. By default,
				all genes are reported when list is empty.
			@DataDir = string of the input directory of where the data is 
    '''
    if 'CCLE' in Dataset:
        df = pd.read_csv(DataDir + Dataset + '/CCLE_expression.csv', sep=',')
        df = pd.melt(df, id_vars=['Unnamed: 0']).rename(columns={'Unnamed: 0': 'Barcode','variable':'GeneName','value':'Value'})
        df['GeneName'] = df['GeneName'].str.split(' ', expand=True)[0]
    elif 'TCGA' in Dataset:
        ListOfFiles = glob.glob(os.path.join(DataDir + Dataset + "/Xena/*/HiSeqV2.gz"))
        output = []
        for FileName in ListOfFiles:
            output.append(pd.read_csv(FileName, sep='\t').set_index('sample'))
        df = pd.melt(pd.concat(output, axis=1).reset_index(), id_vars=['sample']).rename(
			columns={'sample':'GeneName', 'value':'Value','variable':'Barcode'})
    if len(GeneSet) != 0:
        df = df[df['GeneName'].isin(GeneSet)]
    else:
        return(df)


def GetPatientAge(DataDir=os.getcwd() + '/Data/Raw/Clinical/'):
	'''
	Returns a dataframe of two columns: 
		`ShortBarcode` = unique barcode IDs for each patient in TCGA with sample identifiers removed
		`age` = patient age at diagnosis
	Args: @DataDir = string of the input directory of where the data is 
	'''
	age = pd.read_csv(DataDir + 'TCGA/TCGA-CDR-SupplementalTableS1.txt', sep='\t')
	age= age[['bcr_patient_barcode','age_at_initial_pathologic_diagnosis']].rename(columns={
			'bcr_patient_barcode':'ShortBarcode', 'age_at_initial_pathologic_diagnosis':'age'})
	return(age)


def GetDrugResponseData(Screen='primary', AllDrugs=False, DataDir = os.getcwd() + "/Data/Raw/Drug/"):
	'''
	Returns raw cell viability values in the drug screen data (PRISM) from DepMap.
	Args: @Screen = a string for which dataset to read in, `primary` or `secondary` screen in PRISM.
	      @AllDrugs = boolean that when False only outputs viabilities from drugs only targeting proteostasis complexes.
		  @DataDir = string of the input directory of where the data is 
	'''
	drug = pd.read_csv(DataDir + Screen + '-screen-replicate-collapsed-logfold-change.csv', sep=',').rename(
		columns={'Unnamed: 0':'Barcode'}).melt(id_vars='Barcode')  
	drug[['column_name','dose','screenType']] = drug['variable'].str.split('::', expand=True)[[0,1,2]]
	info = pd.read_csv(DataDir + Screen + '-screen-replicate-collapsed-treatment-info.csv', sep=',')
	if AllDrugs: 
		DrugTargets = info[['broad_id','moa']].dropna().rename(columns={'moa':'subgroup'}).assign(Group='All')
	else: # Only look at proteostasis inhibitors (proteasome, hsp90, etc)
		DrugTargets =  pd.concat([
			pd.DataFrame({'broad_id': info[info['moa'] == 'protein synthesis inhibitor']['broad_id'], 
				'Group': 'Translation', 'subgroup':'Protein Synthesis Inhibitor'}),
			pd.DataFrame({ 'broad_id': info[info['moa'] == 'RNA synthesis inhibitor']['broad_id'], 
				'Group': 'Translation', 'subgroup':'RNA Synthesis Inhibitor'}),
			pd.DataFrame({ 'broad_id': info[(info['target'].str.contains('HSP90', na=False)) & (info['moa'].str.contains('HSP inhibitor', na=False))]['broad_id'],  
				'Group': 'Chaperone', 'subgroup':'HSP90 Inhibitor'}),
			pd.DataFrame({ 'broad_id': info[(info['moa'].str.contains('proteasome', na=False))]['broad_id'],
				'Group': 'Proteasome', 'subgroup':'Proteasome Inhibitor'}),
			pd.DataFrame({ 'broad_id': info[(info['moa'] == 'ubiquitin specific protease inhibitor')]['broad_id'],  
				'Group': 'Ubiquitin', 'subgroup':'Ubiquitin-Specific Proteasome Inhibitor'})])
	info = info[['column_name','name']].drop_duplicates()
	drug = info.merge(drug, left_on='column_name', right_on='variable')
	drug = drug.merge(DrugTargets, left_on='column_name_y', right_on='broad_id').rename(columns={'value':'Value'})
	return(drug)


def GetRNAiResponseData(GeneSet = [], DataDir = os.getcwd() + "/Data/Raw/Drug/"):
    '''
	Returns cell viability values from the shRNA screen data (Achilles) from DepMap in the form of a dataframe 
	with 3 columns: Barcode, GeneName, Value.	
	Args: 	@GeneSet = optional argument of a list of gene names to only shRNA values for. By default,
				all genes are reported when list is empty.
			@DataDir = string of the input directory of where the data is 
    '''
    DataDir = os.getcwd() + "/Data/Raw/RNAi/"
    ko = pd.read_csv( DataDir + 'D2_Achilles_gene_dep_scores.csv', sep=',')
    ko = pd.melt(ko, id_vars=['Unnamed: 0']).rename(columns={'Unnamed: 0': 'GeneName','variable':'Barcode','value':'Value'})
    info = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/CCLE/sample_info.csv', sep=',')[['CCLE_Name','DepMap_ID']]
    ko = ko.merge(info, left_on='Barcode',right_on='CCLE_Name').dropna()[['GeneName','Value','DepMap_ID']].rename(columns={'DepMap_ID':'Barcode'})
    ko['GeneName'] = ko['GeneName'].str.split(' ', expand=True)[0]
    if len(GeneSet) != 0:
        ko = ko[ko['GeneName'].isin(GeneSet)]
    return(ko)		



def GetPointMutations(Dataset, DataDir = os.getcwd() + '/Data/Raw/Mutations/'):
	'''
	Returns type and gene name of each mutation in different datasets of a dataframe with 3 columns: 
		Barcode, Variant_Classificaion, Hugo_Symbol.
	Args: 	@Dataset = a string of which dataset to read values from ('CCLE', or 'TCGA')
			@DataDir = string of the input directory of where the data is 
	'''
	if 'TCGA' in Dataset:
		muts = pd.read_csv(DataDir + Dataset + "/mc3.v0.2.8.PUBLIC.maf.gz", sep='\t', 
			usecols = ['Variant_Classification','Tumor_Sample_Barcode','Hugo_Symbol']).rename(columns={'Tumor_Sample_Barcode': 'Barcode'})
		muts['Barcode'] = muts['Barcode'].str[0:15]
	elif 'CCLE' in Dataset:
		muts = pd.read_csv(DataDir + Dataset +  "/CCLE_mutations.csv", sep=',', 
		usecols = ['DepMap_ID','Variant_Classification','Hugo_Symbol']).rename(columns={'DepMap_ID': 'Barcode'})
	return(muts)


def GetCopyNumberMutations(Dataset, DataDir = os.getcwd() + '/Data/Raw/Mutations/'):
	'''
	Returns counts of the number of deletions and amplifications in different datasets in the form of a dataframe with 3 columns: 
		`SAMPLE_NAME` = string identifying each unique barcode ()
		`gain` = counts of the number of amplifications in each sample
		`loss` = counts of the number of deletions in each sample
	Args: 	@Dataset = a string of which dataset to read values from ('CCLE', or 'TCGA')
			@DataDir = string of the input directory of where the data is 
	'''
	if Dataset == 'TCGA':
		SampleInfo = pd.read_csv(DataDir + Dataset + '/CosmicSample.tsv.gz', sep='\t') # Downloaded on 04/07/2020 v91
		cnvs = pd.read_csv(DataDir + Dataset + '/CosmicCompleteCNA.tsv.gz', sep='\t', usecols=['SAMPLE_NAME','MUT_TYPE'])
		cnvs = cnvs[cnvs['SAMPLE_NAME'].str.contains('TCGA', na=False)]
	elif Dataset == 'CCLE':
		cnvs = pd.read_csv(DataDir + Dataset + '/CCLE_ABSOLUTE_combined_20181227.txt', sep='\t', # Version CCLE 2019 in DepMap
				usecols=['depMapID','Modal_HSCN_1','Modal_HSCN_2','Modal_Total_CN', 'LOH','Homozygous_deletion']).rename(
					columns={'depMapID': 'SAMPLE_NAME'}).drop_duplicates()
		cnvs['MUT_TYPE'] = np.nan
		cnvs = cnvs[~((cnvs['Modal_Total_CN'] == 2) & (cnvs['Modal_HSCN_1'] == 1) & (cnvs['Modal_HSCN_2'] == 1))] ### remove probes that are a CN of 2
		cnvs.loc[cnvs['Modal_Total_CN'] > 2, 'MUT_TYPE'] = 'gain' ### annotate amplification similar to COSMIC notation
		cnvs.loc[cnvs['Modal_Total_CN'] < 2, 'MUT_TYPE'] = 'loss' ### annotate deletion similar to COSMIC notation
	return(cnvs.groupby(['SAMPLE_NAME'])['MUT_TYPE'].value_counts().reset_index(name='counts').pivot(
			index='SAMPLE_NAME', columns='MUT_TYPE', values='counts').reset_index().replace(np.nan,0))


def AnnotateMutationalLoad(muts, MutType):
	'''
	Takes in raw mutation calls and calculates number of mutations of a mutation type of interest in each sample.
	Args:	@muts = a dataframe of raw mutation calls (SNVs or CNVs) from either @GetCopyNumberMutations or @GetPointMutations
			@MutType = a string of which mutation calls to filter/aggregate for 
	'''
	if MutType == 'KsKa': # just protein coding mutations
		MutationClasses = dict(Ka={'Missense_Mutation', 'Nonsense_Mutation'},  Ks={'Silent'})
		KaOrKs = {k:muts['Variant_Classification'].isin(frozenset(v)) for k, v in MutationClasses.items()}
		muts = pd.concat([muts, pd.DataFrame(KaOrKs) * 1], axis=1)
		muts['TotalKaKs'] = muts['Ka'] + muts['Ks']
		MutationRate = muts.groupby('Barcode')['TotalKaKs'].sum().reset_index().rename(columns={'TotalKaKs':'MutLoad'})
	elif MutType == 'CNV': # total number of amplifications and snvs
		MutationRate = muts.set_index(['SAMPLE_NAME']).sum(axis=1).reset_index().rename(columns={0:'MutLoad', 'SAMPLE_NAME':'Barcode'})
	elif MutType == 'CN_Deletions':
		MutationRate = muts[['SAMPLE_NAME','gain']].rename(columns={'SAMPLE_NAME':'Barcode','gain':'MutLoad'})
	elif MutType == 'CN_Amplifications':
		MutationRate = muts[['SAMPLE_NAME','loss']].rename(columns={'SAMPLE_NAME':'Barcode','loss':'MutLoad'})
	elif MutType == 'SNV': 
		muts = muts[~muts['Variant_Classification'].isin(['In_Frame_Ins','In_Frame_Del'])] #remove MNVs
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'Nonsynonymous':
		muts = muts[muts['Variant_Classification'] == 'Missense_Mutation']
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'Synonymous':
		muts = muts[muts['Variant_Classification'] == 'Silent']
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'Nonsense':
		muts = muts[muts['Variant_Classification'] == 'Nonsense_Mutation']
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	return(MutationRate)


def GetTumorPurity(DataDir = os.getcwd() + '/Data/Raw/Mutations/TCGA/' ):
    '''Returns a dataframe of tumor purity calls for each patient in TCGA with the following columns: `Barcode`,`purity`'''
    purity = pd.read_csv(DataDir + 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t')
    purity['Tumor_Sample_Barcode'] = purity['array']
    purity = purity[['Tumor_Sample_Barcode', 'purity']].rename(columns={'Tumor_Sample_Barcode':'Barcode'})
    purity['Barcode'] = purity['Barcode'].str[0:15]
    return(purity)


def GetTissueType(Dataset, DataDir = os.getcwd() + '/Data/Raw/Mutations/' ):
	'''
	Returns a dataframe of the cancer type associated with each tumor/cell line with the following columns: `Barcode`,`type`
	Args: 	@Dataset = a string of which dataset to read expression values from ('CCLE', or 'TCGA')
			@DataDir = string of the input directory of where the data is 
	'''
	if 'TCGA' in Dataset:
		muts = pd.read_csv(DataDir + Dataset + '/CancerTypesTCGA', sep='\t').rename(
			columns={'Tumor_Sample_Barcode': 'Barcode','subtype':'type'})[['Barcode','type']]		
		muts['Barcode'] = muts['Barcode'].str[0:15]
	elif 'CCLE' in Dataset:
		muts = pd.read_csv(DataDir +  Dataset + '/sample_info.csv', sep=',', usecols=[ 'disease', 'DepMap_ID']).rename(
			columns={'DepMap_ID': 'Barcode','disease':'type'})
	return(muts)



def GetTopMutationalLoadTumors(TopPercentThreshold = 0.1, PerCancerType=True):
	df = AnnotateMutationalLoad(GetPointMutations('TCGA'), 'KsKa').merge(GetTissueType('TCGA'), left_on='Barcode', right_on='Barcode')
	if PerCancerType:
		topload = df.groupby('type').apply(lambda x: x.nlargest(int(len(x) * TopPercentThreshold), 'MutLoad')).reset_index(drop=True)
		percent = topload.groupby('type')['MutLoad'].min().reset_index().assign(MutationLoadPercentile = str(int((TopPercentThreshold * 100))) + '%')
	else: # Top load by all cancer types
		topload = df['MutLoad'].nlargest(int(len(df) * TopPercentThreshold)).reset_index().drop(columns='index')
		percent = pd.DataFrame({'MutLoad':topload['MutLoad'].min()}, index=[0]).assign(
			MutationLoadPercentile = str(int((TopPercentThreshold * 100))) + '%').assign(type='All Cancer Types')
	return(percent.rename(columns={'type':'Cancer Type','MutLoad':'Minimum Threshold of Mutations', 'MutationLoadPercentile': 'Mutation Load Percentile'}))
