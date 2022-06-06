'''
This script contains gene annotations used to determine proteostasis gene groups and
drug annotations used to group drugs into broad functional categories.
'''

import pandas as pd

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