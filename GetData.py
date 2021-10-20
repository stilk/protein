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
        df = pd.read_csv('TCGA_TS.txt')
    if len(GeneSet) != 0:
        df = df[df['GeneName'].isin(GeneSet)]
    else:
        return(df)


def GetDrugResponseData(Screen):
    '''
    Reads in drug screen data (PRISM) from DepMap.
    Args: @screen = a string for which dataset to read in, `primary` or `secondary` screen
    '''
    DataDir = os.getcwd() + "/Data/Raw/Drug/"
    drug = pd.read_csv(DataDir + Screen + '-screen-replicate-collapsed-logfold-change.csv', sep=',').rename(
        columns={'Unnamed: 0':'Barcode'}).melt(id_vars='Barcode')  
    drug[['column_name','dose','screenType']] = drug['variable'].str.split('::', expand=True)[[0,1,2]]
    info = pd.read_csv(DataDir + Screen + '-screen-replicate-collapsed-treatment-info.csv', sep=',')
    info = info[['column_name','name']].drop_duplicates()
    drug = info.merge(drug, left_on='column_name', right_on='variable')
    #return(drug[druge['column_name_y'].isin(pd.Series(drugNames).str.split('::', expand=True)[0])])
    return(drug)


def GetRNAiResponseData(GeneSet = []):
    DataDir = os.getcwd() + "/Data/Raw/RNAi/"
    ko = pd.read_csv( DataDir + 'D2_combined_gene_dep_scores.csv', sep=',')
    ko = pd.melt(ko, id_vars=['Unnamed: 0']).rename(columns={'Unnamed: 0': 'GeneName','variable':'Barcode','value':'Value'})
    info = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/CCLE/sample_info.csv', sep=',')[['CCLE_Name','DepMap_ID']]
    ko = ko.merge(info, left_on='Barcode',right_on='CCLE_Name').dropna()[['GeneName','Value','DepMap_ID']].rename(columns={'DepMap_ID':'Barcode'})
    ko['GeneName'] = ko['GeneName'].str.split(' ', expand=True)[0]
    if len(GeneSet) != 0:
        ko = ko[ko['GeneName'].isin(GeneSet)]
    return(ko)		


def DownloadFastaForVep():
	'''
	To run the --hgvs tag with VEP, a fasta file is required. This function downloads the human
	genome (GRCh37) Fasta file from Ensembl. For documentation see:
	https://uswest.ensembl.org/info/docs/webcode/mirror/tools/vep.html
	'''
	VEPDir='/labs/ccurtis2/tilk/software/vep/'
	if os.path.isfile(VEPDir + 'Homo_sapiens.GRCh37.75.dna.toplevel.bgzip.fa.gz') == False:
		subprocess.call(['bash', '-c', 
		'wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz' +
		'-P ' + VEPDir ])
		# Fasta must be compressed with bgzip, not gzip
		subprocess.call(['bash', '-c', 'zcat ' + VEPDir + 'Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz | bgzip -c > ' +
			VEPDir + 'Homo_sapiens.GRCh37.75.dna.toplevel.bgzip.fa.gz'])
	else:
		return(print('Fasta already downloaded.'))

def DownloadVEPIndex():
	'''
	VEP requires a cached index to run mutation annotations. This function downloads the GRCh37 genom
	index to your home directory. This downloads a large file (>10G) and takes a very long time run (~9 hrs). 
	'''	
	import os
	import subprocess
	if os.path.isfile('~/.vep/homo_sapiens_vep_104_GRCh37.tar.gz') == False: # If VEP output doesn't exist yet
		# Download cached index
		subprocess.call(['bash', '-c', 
		'wget http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/homo_sapiens_vep_104_GRCh37.tar.gz -P /labs/ccurtis2/tilk/software/vep'])
		# Un-tar package
		subprocess.call(['bash','-c', 'tar -zxvf /labs/ccurtis2/tilk/software/vep/homo_sapiens_vep_104_GRCh37.tar.gz'])
		return(print('Finished downloading index.'))
	else:
		return(print('Index already created.'))


def RunVEPToGetDeleteriousScores():
	'''
	Runs VEP (Variant Effect Predictor) to obtain polyphen/sift scores for CCLE
	Required columns to run VEP are: Chromosome, Start, Stop, Ref/Alt, Strand
	'''
	muts = pd.read_csv( os.getcwd() + "/Data/Raw/Mutations/CCLE/CCLE_mutations.csv")
	# Generate VEP Input to file
	muts['Ref/Alt'] = muts['Reference_Allele'] + '/' + muts['Tumor_Seq_Allele1']
	VEPInput = muts[['Chromosome','Start_position','End_position','Ref/Alt','Strand']]
	VEPInput.to_csv(os.getcwd() + "/Data/Raw/Mutations/CCLE/CCLE_VEP_Input" , sep= ' ', index=False, header=None)
	# Run VEP
	if os.path.isfile(os.getcwd() + "/Data/Raw/Mutations/CCLE/CCLE_VEP_Output") == False: # If VEP output doesn't exist yet
		#DownloadVEPIndex()
		subprocess.call(['bash', '-c', 
			'vep -e -v --database --dir /labs/ccurtis2/tilk/software/vep/ --species homo_sapiens -i ' + 
			os.getcwd() + "/Data/Raw/Mutations/CCLE/CCLE_VEP_Input" + ' -o ' + os.getcwd() + "/Data/Raw/Mutations/CCLE/CCLE_VEP_Output" +
			'--sift --polyphen --format guess --force_overwrite'])
		VEPOutput = pd.read_csv(os.getcwd() + "/Data/Raw/Mutations/CCLE/CCLE_VEP_Output")
	return(pd.concat([muts, VEPOutput], axis=1))

def GetStability(Dataset):
	DataDir = "/labs/ccurtis2/tilk/scripts/protein/Data/MutFunc/"
	if Dataset == 'TCGA':
		return(pd.read_csv(DataDir + 'TCGA_ProteinStability_HomologyModel'))

def GetMutations(Dataset):
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
    elif 'GTEX' in Dataset:
        muts = pd.read_csv(DataDir + 'Table_mutations_all_Garcia_Nieto_et_al_Genome_Biology.tsv', sep='\t', usecols = ['gtexIds_samples']).rename(
            columns={'gtexIds_samples':'Barcode'})
    return(muts)



def AnnotateMutationalLoad(muts, MutType):
	if MutType == 'KsKa': # just protein coding mutations
		MutationClasses = dict(Ka={'Missense_Mutation', 'Nonsense_Mutation'},  Ks={'Silent'})
		KaOrKs = {k:muts['Variant_Classification'].isin(frozenset(v)) for k, v in MutationClasses.items()}
		muts = pd.concat([muts, pd.DataFrame(KaOrKs) * 1], axis=1)
		muts['TotalKaKs'] = muts['Ka'] + muts['Ks']
		MutationRate = muts.groupby('Barcode')['TotalKaKs'].sum().reset_index().rename(columns={'TotalKaKs':'MutLoad'})
	elif MutType == 'CNV_Altered':
		chrom = pd.read_csv('/labs/ccurtis2/tilk/02_DNDS/CN_Analysis/02_geneSets/annotations_dEdI/biomart_chromosomes.csv' , sep=',')
		MutationRate = pd.DataFrame(cn.groupby('Barcode')['Length'].sum()/(chrom['Length'].sum().astype(int))).reset_index().rename(columns={'Length':'MutLoad'})
	elif MutType == 'Deleterious':
		MutationRate = muts[muts['Deleterious_Meta'] == True].groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'SNV': # all SNV mutations
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'Nonsynonymous':
		muts = muts[muts['Variant_Classification'] == 'Missense_Mutation']
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'Polyphen':
		muts['PolyPhen'] = muts['PolyPhen'].str.split("(", expand=True)[0]
		muts = muts[muts['PolyPhen'].isin(['probably_damaging', 'possibly_damaging'])]
		MutationRate = muts.groupby('Barcode').size().reset_index().rename(columns={0:'MutLoad'})
	elif MutType == 'NonsynonymousHE':
		print('foo')
	return(MutationRate)

def GetHighlyExpressedGenes():
    ''' Defined as having TPM greater than 1000 across all samples in GTEX.'''
    OutDir='/labs/ccurtis2/tilk/scripts/cancer-HRI/Data/GeneSets/'
    samp=pd.read_csv(OutDir + 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', skiprows=2, sep='\t')
    samp['Name'] = samp['Name'].str.split('.', expand=True)[0]
    genes = samp.set_index(['Description','Name'])
    return(genes.loc[(genes > 500).all(axis=1)].reset_index()[['Description','Name']])


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
		muts = GetMutations(Dataset='TCGA')
		mut_rate = AnnotateMutationalLoad(muts, MutType='SNV')
		tissue = GetTissueType(Dataset)[['Barcode','type']]
		tissue.merge(mut_rate, left_on='Barcode' ,right_on='Barcode').to_csv(
			os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE'
		)

