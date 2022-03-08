import pandas as pd
import numpy as np
import glob
import os
from scipy import stats
import scipy
import subprocess

def GetExpressionData(Dataset, GeneSet=[]):
    '''
    Returns a data frame with 3 columns: Barcode, GeneName, Value
    Optional argument of @GeneSet, where gene expression values are subsetted for a list of genes.
    '''
    DataDir = os.getcwd() + "/Data/Raw/Expression/" + Dataset
    if 'CCLE' in Dataset:
        df = pd.read_csv(DataDir + '/CCLE_expression.csv', sep=',')
        df = pd.melt(df, id_vars=['Unnamed: 0']).rename(columns={'Unnamed: 0': 'Barcode','variable':'GeneName','value':'Value'})
        df['GeneName'] = df['GeneName'].str.split(' ', expand=True)[0]
    elif 'GTEX' in Dataset:
        df = pd.read_csv(DataDir + '/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2).drop(columns=['Name'])
        df = pd.melt(df, id_vars=['Description']).rename(columns={'Description': 'GeneName','variable':'Barcode','value':'Value'})
    elif 'TCGA' in Dataset:
        ListOfFiles = glob.glob(os.path.join(DataDir + "/Xena/*/HiSeqV2.gz"))
        output = []
        for FileName in ListOfFiles:
            output.append(pd.read_csv(FileName, sep='\t').set_index('sample'))
        df = pd.melt(pd.concat(output, axis=1).reset_index(), id_vars=['sample']).rename(columns={'sample':'GeneName', 'value':'Value','variable':'Barcode'})
    elif 'TS' in Dataset:
        df = pd.read_csv('/oak/stanford/scg/lab_ccurtis2/tilk/scripts/protein/TCGA_TS.txt')
    if len(GeneSet) != 0:
        df = df[df['GeneName'].isin(GeneSet)]
    else:
        return(df)



def GetPatientAge(Dataset='TCGA'):
	DataDir=os.getcwd() + '/Data/Raw/Clinical/'
	if Dataset == 'TCGA':
		age = pd.read_csv(DataDir + 'TCGA/TCGA-CDR-SupplementalTableS1.txt', sep='\t')
		age= age[['bcr_patient_barcode','age_at_initial_pathologic_diagnosis']].rename(columns={
			'bcr_patient_barcode':'ShortBarcode', 'age_at_initial_pathologic_diagnosis':'age'})
	return(age)


def GetProteinExpressionData(Dataset, GeneSet=[]):
	'''
	Returns a data frame columns with (mainly) 3 columns: Barcode, and value, geneName.
	Optional argument of geneSet, where you get protein expression values only for a list of genes.
	'''
	if Dataset == 'CCLE':
		protein = pd.read_csv('/labs/ccurtis2/tilk/07_PRISM/00_RawData/Protein/protein_quant_current_normalized.csv.gz', sep=',')
		protein = protein.drop(['Protein_Id','Description','Group_ID','Uniprot','Uniprot_Acc'], axis=1).set_index('Gene_Symbol').transpose().unstack().reset_index().rename(columns={'level_1':'Protein_Id',0:'ProteinQuant'})
		protein[['Protein_Id','Replicate']] = protein.Protein_Id.str.rsplit('_', 1, expand=True) ### split into replicates and IDs
		protein = protein[protein['Replicate'] != 'Peptides'] ### remove weird samples for now 
		info = pd.read_csv('/labs/ccurtis2/tilk/07_PRISM/00_RawData/SNV/sample_info.csv', sep=',')[['CCLE_Name','DepMap_ID']]
		df = protein.merge(info, left_on='Protein_Id',right_on='CCLE_Name').dropna().rename(columns={'DepMap_ID':'Barcode','Gene_Symbol': 'GeneName','ProteinQuant':'Value'})
		if len(GeneSet) != 0:
				df = df[df['GeneName'].isin(GeneSet)]
		return(df)


def GetDrugResponseData(Screen='primary', AllDrugs=False):
	'''
	Reads in drug screen data (PRISM) from DepMap.
	Args: @screen = a string for which dataset to read in, `primary` or `secondary` screen
	      @AllDrugs = boolean of whether to only output protostatis drugs
	'''
	DataDir = os.getcwd() + "/Data/Raw/Drug/"
	drug = pd.read_csv(DataDir + Screen + '-screen-replicate-collapsed-logfold-change.csv', sep=',').rename(
		columns={'Unnamed: 0':'Barcode'}).melt(id_vars='Barcode')  
	drug[['column_name','dose','screenType']] = drug['variable'].str.split('::', expand=True)[[0,1,2]]
	info = pd.read_csv(DataDir + Screen + '-screen-replicate-collapsed-treatment-info.csv', sep=',')
	if AllDrugs: 
		DrugTargets = info[['broad_id','moa']].dropna().rename(columns={'moa':'subgroup'}).assign(Group='All')
	else: # Only look at proteostasis inhibitors (proteasome, hsp90, etc)
		DrugTargets =  pd.concat([pd.DataFrame({'broad_id': info[info['moa'] == 'protein synthesis inhibitor']['broad_id'], 'Group': 'Translation', 'subgroup':'Protein Synthesis Inhibitor'}),
								pd.DataFrame({'broad_id': info[info['moa'] == 'RNA synthesis inhibitor']['broad_id'], 'Group': 'Translation', 'subgroup':'RNA Synthesis Inhibitor'}),
								pd.DataFrame({ 'broad_id': info[(info['target'].str.contains('HSP90', na=False)) & (info['moa'].str.contains('HSP inhibitor', na=False))]['broad_id'],  'Group': 'Chaperone', 'subgroup':'HSP90 Inhibitor'}),
								pd.DataFrame({ 'broad_id': info[(info['moa'].str.contains('proteasome', na=False))]['broad_id'],  'Group': 'Proteasome', 'subgroup':'Proteasome Inhibitor'}),
								pd.DataFrame({ 'broad_id': info[(info['moa'] == 'ubiquitin specific protease inhibitor')]['broad_id'],  'Group': 'Ubiquitin', 'subgroup':'Ubiquitin-Specific Proteasome Inhibitor'})
		])
	info = info[['column_name','name']].drop_duplicates()
	drug = info.merge(drug, left_on='column_name', right_on='variable')
	drug = drug.merge(DrugTargets, left_on='column_name_y', right_on='broad_id').rename(columns={'value':'Value'})
	#drugInfo = pd.read_csv('/labs/ccurtis2/tilk/07_PRISM/00_RawData/Viability/DrugScreen/PrimaryScreen/primary-screen-replicate-collapsed-treatment-info.csv', sep=',')
	#drugTargets = drugInfo.set_index(['broad_id','name'])['target'].str.split(', ', expand=True).reset_index().melt(id_vars=['broad_id', 'name']).dropna()
	#return(drug[druge['column_name_y'].isin(pd.Series(drugNames).str.split('::', expand=True)[0])])
	return(drug)


def GetRNAiResponseData(GeneSet = []):
    '''
    Not all cancer cell lines were profiled evenly in the shRNA screen. Only reports RNAi data for
    genes that have at least half of the median number of barcodes cell lines profiled. 
    '''
    DataDir = os.getcwd() + "/Data/Raw/RNAi/"
    ko = pd.read_csv( DataDir + 'D2_Achilles_gene_dep_scores.csv', sep=',')
	#ko = pd.read_csv( DataDir + 'D2_DRIVE_gene_dep_scores.csv', sep=',')
    ko = pd.melt(ko, id_vars=['Unnamed: 0']).rename(columns={'Unnamed: 0': 'GeneName','variable':'Barcode','value':'Value'})
    info = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/CCLE/sample_info.csv', sep=',')[['CCLE_Name','DepMap_ID']]
    ko = ko.merge(info, left_on='Barcode',right_on='CCLE_Name').dropna()[['GeneName','Value','DepMap_ID']].rename(columns={'DepMap_ID':'Barcode'})
    ko['GeneName'] = ko['GeneName'].str.split(' ', expand=True)[0]
    #NumBarcodesScreened = ko.groupby(['GeneName']).size().reset_index().rename(columns={0:'NumBarcodesScreened'})
    #ko = ko.merge(NumBarcodesScreened, left_on='GeneName', right_on='GeneName')
    #ko = ko[ko['NumBarcodesScreened'] > ko['NumBarcodesScreened'].median()/2] 
    if len(GeneSet) != 0:
        ko = ko[ko['GeneName'].isin(GeneSet)]
    return(ko)		



def GetPointMutations(Dataset):
    DataDir = os.getcwd() + '/Data/Raw/Mutations/' + Dataset + '/'
    if 'TCGA' in Dataset:
        muts = pd.read_csv(DataDir + "mc3.v0.2.8.PUBLIC.maf.gz", sep='\t', usecols = ['Variant_Classification','Tumor_Sample_Barcode','Hugo_Symbol','PolyPhen']).rename(
            columns={'Tumor_Sample_Barcode': 'Barcode'})
        muts['Barcode'] = muts['Barcode'].str[0:15]
        muts = muts[~muts['Variant_Classification'].isin(['In_Frame_Ins','In_Frame_Del'])] # Don't look at MNVs
    elif 'CCLE' in Dataset:
        muts = pd.read_csv(DataDir + "CCLE_mutations.csv", sep=',', usecols = ['DepMap_ID','Variant_Classification','Hugo_Symbol']).rename(
            columns={'DepMap_ID': 'Barcode'})
        muts = muts[~muts['Variant_Classification'].isin(['In_Frame_Ins','In_Frame_Del'])] # Don't look at MNVs
    return(muts)


def GetCopyNumberMutations(Dataset):
	DataDir = os.getcwd() + '/Data/Raw/Mutations/' + Dataset + '/'
	if Dataset == 'TCGA':
		SampleInfo = pd.read_csv(DataDir + 'CosmicSample.tsv.gz', sep='\t') # Downloaded on 04/07/2020 v91
		cnvs = pd.read_csv(DataDir + 'CosmicCompleteCNA.tsv.gz', sep='\t', usecols=['SAMPLE_NAME','MUT_TYPE'])
		cnvs = cnvs[cnvs['SAMPLE_NAME'].str.contains('TCGA', na=False)]
	elif Dataset == 'CCLE':
		cnvs = pd.read_csv(DataDir + 'CCLE_ABSOLUTE_combined_20181227.txt', sep='\t', # Version CCLE 2019 in DepMap
				usecols=['depMapID','Modal_HSCN_1','Modal_HSCN_2','Modal_Total_CN', 'LOH','Homozygous_deletion']).rename(columns={'depMapID': 'SAMPLE_NAME'}).drop_duplicates()
		cnvs['MUT_TYPE'] = np.nan
		cnvs = cnvs[~((cnvs['Modal_Total_CN'] == 2) & (cnvs['Modal_HSCN_1'] == 1) & (cnvs['Modal_HSCN_2'] == 1))] ### remove probes that are a CN of 2
		cnvs.loc[cnvs['Modal_Total_CN'] > 2, 'MUT_TYPE'] = 'gain' ### annotate amplification similar to COSMIC notation
		cnvs.loc[cnvs['Modal_Total_CN'] < 2, 'MUT_TYPE'] = 'loss' ### annotate deletion similar to COSMIC notation
		#cnvs.loc[cnvs['Homozygous_deletion'] == 1, 'MUT_TYPE'] = 'loss' ### annotate homozygous deletion similar to COSMIC notation
	return(cnvs.groupby(['SAMPLE_NAME'])['MUT_TYPE'].value_counts().reset_index(name='counts').pivot(
			index='SAMPLE_NAME', columns='MUT_TYPE', values='counts').reset_index().replace(np.nan,0))




def AnnotateMutationalLoad(muts, MutType):
	if MutType == 'KsKa': # just protein coding mutations
		MutationClasses = dict(Ka={'Missense_Mutation', 'Nonsense_Mutation'},  Ks={'Silent'})
		KaOrKs = {k:muts['Variant_Classification'].isin(frozenset(v)) for k, v in MutationClasses.items()}
		muts = pd.concat([muts, pd.DataFrame(KaOrKs) * 1], axis=1)
		muts['TotalKaKs'] = muts['Ka'] + muts['Ks']
		MutationRate = muts.groupby('Barcode')['TotalKaKs'].sum().reset_index().rename(columns={'TotalKaKs':'MutLoad'})
	elif MutType == 'CNV':
		MutationRate = muts.set_index(['SAMPLE_NAME']).sum(axis=1).reset_index().rename(columns={0:'MutLoad', 'SAMPLE_NAME':'Barcode'})
	elif MutType == 'CN_Deletions':
		MutationRate = muts[['SAMPLE_NAME','gain']].rename(columns={'SAMPLE_NAME':'Barcode','gain':'MutLoad'})
	elif MutType == 'CN_Amplifications':
		MutationRate = muts[['SAMPLE_NAME','loss']].rename(columns={'SAMPLE_NAME':'Barcode','loss':'MutLoad'})
	elif MutType == 'Deleterious':
		MutationRate = muts[muts['Deleterious_Meta'] == True].groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'SNV': # all SNV mutations
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
	elif MutType == 'Polyphen':
		muts['PolyPhen'] = muts['PolyPhen'].str.split("(", expand=True)[0]
		muts = muts[muts['PolyPhen'].isin(['probably_damaging', 'possibly_damaging'])]
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	return(MutationRate)


def GetTumorPurity():
    ''' Reads in a dataframe of mutations and adds a column of tumor purity calls for each patient.'''
    purity = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t')
    purity['Tumor_Sample_Barcode'] = purity['array']
    purity = purity[['Tumor_Sample_Barcode', 'purity']].rename(columns={'Tumor_Sample_Barcode':'Barcode'})
    purity['Barcode'] = purity['Barcode'].str[0:15]
    return(purity)


def GetTissueType(Dataset):
	DataDir = os.getcwd() + '/Data/Raw/Mutations/' + Dataset + '/'
	if 'TCGA' in Dataset:
		muts = pd.read_csv(DataDir + "CancerTypesTCGA", sep='\t').rename(columns={'Tumor_Sample_Barcode': 'Barcode','subtype':'type'})[['Barcode','type']]		
		muts['Barcode'] = muts['Barcode'].str[0:15]
	elif 'GTEX' in Dataset:
		muts = pd.read_csv(DataDir + "Table_mutations_all_Garcia_Nieto_et_al_Genome_Biology.tsv", sep='\t', 
		    usecols=['gtexIds_samples','tissue']).rename(columns={'gtexIds_samples': 'Barcode','tissue':'type'}).drop_duplicates()
	elif 'CCLE' in Dataset:
		muts = pd.read_csv(DataDir + 'sample_info.csv', sep=',', usecols=[ 'disease', 'DepMap_ID']).rename(columns={'DepMap_ID': 'Barcode','disease':'type'})
	return(muts)

    
def GetNumberOfMutsAndCancerType(Dataset='TCGA'):
	if Dataset == 'TCGA':
		muts = GetPointMutations(Dataset='TCGA')
		mut_rate = AnnotateMutationalLoad(muts, MutType='KsKa')
		tissue = GetTissueType(Dataset)[['Barcode','type']]
		tissue.merge(mut_rate, left_on='Barcode' ,right_on='Barcode').to_csv(
			os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE'
		)



