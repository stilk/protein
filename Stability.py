

import pandas as pd
import numpy as np
import glob
import os

def GetProteinStabilityInfo(Model):
    '''
    Stability calls were downloaded from MutFunc.
    Uniprot IDs were downloaded from: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

    '''
    uniprot = pd.read_csv(os.getcwd() + '/Data/MutFunc/HUMAN_9606_idmapping.dat.gz', sep='\t', header=None)
    uniprot = uniprot[uniprot[1].isin(['Ensembl'])].rename(columns={0:'Uniprot', 2:'Gene'}) # Link Ensembl gene IDs (ENSG) to Uniprot
    if Model == 'Experimental_Model':
        df = pd.read_csv('/labs/ccurtis2/tilk/09_PROTEIN/08_Stability/exp.tab.gz', sep='\t')
    elif Model == 'Homology_Model':
        df = pd.read_csv('/labs/ccurtis2/tilk/09_PROTEIN/08_Stability/mod.tab.gz', sep='\t')
    df = df.merge(uniprot[['Uniprot','Gene']].drop_duplicates(), left_on='uniprot_id', right_on='Uniprot') ### adds ENSG IDs to stability info
    df['uniprot_pos'] = df['uniprot_pos'].astype(str)
    return(df)


def MatchStabilityToTCGA(Model):
    muts = pd.read_csv("/labs/ccurtis2/tilk/02_DNDS/mutationAnalysis/00_mafFiles/mc3/mc3.v0.2.8.PUBLIC.maf.gz", sep='\t')
    muts['MutID'] = list(range(0, len(muts)))
    muts[['aa_wt','uniprot_pos','aa_mt']] = muts['HGVSp_Short'].str.split('p.', expand=True)[1].str.extract('([A-Za-z]+)(\d+\.?\d*)([A-Za-z]*)', expand = True)
    df = GetProteinStabilityInfo(Model)
    stability = muts[['aa_wt','uniprot_pos','aa_mt','Feature','MutID','Tumor_Sample_Barcode','Hugo_Symbol','Gene','Variant_Classification']].rename(
        columns={'Feature':'ENST'}).merge(df, left_on= ['aa_wt','uniprot_pos','aa_mt','Gene'], right_on=['aa_wt','uniprot_pos','aa_mt','Gene'],  how='left')
    return(stability)
    




MatchStabilityToTCGA(Model='Homology_Model').to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/MutFunc/TCGA_ProteinStability_HomologyModel')

