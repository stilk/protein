'''
This script performs all alternative splicing analysis used in the manuscript. 
'''

from GetData import *

def GetAlternativeSplicingData(AS='RI', BarcodesOfInt=[]):
	'''
	Aggregates raw PSI counts for each gene and patient from TCGASpliceSeq by alternative splicing event.
	Args:   @AS = a string for which alternative splicing event to look at (e.g. 'RI' for intron retention)
	        @BarcodesOfInt = a vector of barcodes to only aggregate AS events for. Used in the delta PSI comparison
	            between low and high mutational load tumors.
	'''
	DataDir = '/labs/ccurtis2/tilk/09_PROTEIN/splicingAnalysis/spliceSeq/newDownload/'
	ListOfFiles = glob.glob(os.path.join(DataDir, "PSI_download_*zip"))
	OutDF = pd.DataFrame()
	DataColumns = ['symbol', 'as_id', 'splice_type', 'exons', 'from_exon', 'to_exon']
	# Read and merge alternative splicing calls for each cancer type
	for FileName in ListOfFiles:
		CancerType = pd.read_csv(FileName, sep='\t')
		CancerType = CancerType[CancerType['splice_type'] == AS] # Select AS call types (i.e intron retention, etc)
		if len(BarcodesOfInt) == 0: # Don't subset, choose all barcodes in data frame
			OverlappingBarcodes = CancerType.columns[CancerType.columns.str.contains('TCGA')].tolist()
		else: # Choose barcodes of interest that overlap data frame
			OverlappingBarcodes = pd.Series(CancerType.columns)[pd.Series(CancerType.columns).isin(BarcodesOfInt)].tolist()
			if len(OverlappingBarcodes) == 0: # No barcodes match for this cancer type
				continue
		CancerType = CancerType[OverlappingBarcodes + DataColumns]
		if (len(OutDF) == 0):
			OutDF = CancerType
		else:
			OutDF = OutDF.merge(CancerType, left_on=DataColumns, right_on=DataColumns, how='outer')
		print('Done analyzing: ' + FileName)
	OutDF = OutDF.set_index(DataColumns)
	return(OutDF)


def GetAvgTranscriptLevelsForAS():
    '''
    Outputs a dataframe for each patient whether the expression of each gene is under-expressed or over-expressed relative to 
    that patient's cancer subtype. Genes are counted as ender or over-expressed if they are 1 STD from the mean.
    Output of table is used in @GetExpLevelsForASEvents.
    '''
    Exp = GetExpressionData(Dataset='TCGA')
    Tissue = GetTissueType(Dataset='TCGA')
    Exp = Exp.merge(Tissue, left_on='Barcode', right_on='Barcode')
    Avg = Exp.groupby(['GeneName','type'])['Value'].agg(['mean','std']).rename(columns={'mean':'AvgValue'}).reset_index()# Avg gene expression for each cancer type
    Exp = Exp.merge(Avg, left_on=['GeneName','type'], right_on=['GeneName','type'])
    Exp['DeviationFromMean'] = Exp['Value'] - Exp['AvgValue']
    Exp['DepletedTranscript'] = (Exp['DeviationFromMean'] < 0) & (Exp['DeviationFromMean'].abs() > Exp['std'])
    Exp['OverExpressedTranscript'] = (Exp['DeviationFromMean'] > 0) & (Exp['DeviationFromMean'].abs() > Exp['std'])
    Exp.to_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/' + Dataset + '_AvgGeneExpForAllCancerTypes')
    return(Exp)


def GetExpLevelsForASEvents(ASType, PSI_Threshold, FilterForeQTLS=True):
    '''
    Compares expression levels for each gene by cancer subtype in TCGA and merges alternative splicing calls
    for the same gene. Outputs a table for each patient (1) how many total transcripts there are with AS calls in the
    `True` column and (2) how many of those transcripts are under-expressed in `False` column. To be counted in the `True`
    column, the patient can not have a point mutation in the same gene.
    Parameters:
    @ASType = a string for which alternative splicing event to look at (e.g. 'RI' for intron retention)
    @PSI_Threshold = a float of the threshold of PSI values to be used to determine whether there is an alternative
                    splicing event in a gene. Higher PSI values have more transcripts that contain the AS event.
    '''
    OutDir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    Exp = pd.read_csv('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_AvgGeneExpForAllCancerTypes') #Generated by GetAvgTranscriptLevelsForAS(Dataset='TCGA')
    # Merge expression values for each patient and whether they have a mutation in same AS event
    Exp = Exp.merge(GetPointMutations('TCGA'), left_on=['Barcode','GeneName'], right_on=['Barcode','Hugo_Symbol'], how='left')
    Exp['Hugo_Symbol'] = Exp['Hugo_Symbol'].replace(np.nan, 0) 
    Exp['Barcode'] = Exp['Barcode'].str.replace('-','_').str[0:12]
    # Merge expression values with AS events and use PSI threshold to count events as spliced in or not
    AS = GetAlternativeSplicingData(ASType).stack().reset_index().rename(columns={0:'PSI','level_6':'Barcode','symbol':'GeneName'})
    Out = AS.merge(Exp[['Barcode','GeneName','DeviationFromMean','DepletedTranscript','OverExpressedTranscript','Hugo_Symbol']],
        left_on=['Barcode','GeneName'], right_on=['Barcode','GeneName'])
    Out = pd.concat([Out[Out['PSI'] < PSI_Threshold].assign(Group='Not_Spliced_In'), Out[Out['PSI'] >= PSI_Threshold].assign(Group='Spliced_In')])
    # Save output for supplemental figure to asses how many events are being removed due to potential eQTL effects
    Out[Out['Hugo_Symbol'] != 0].to_csv(OutDir + 'TCGA_AS_EventsRemoved'+ ASType +'_Counts_ThresholdByPSI_' + str(PSI_Threshold))
    Out[Out['Hugo_Symbol'] == 0].to_csv(OutDir + 'TCGA_AS_EventsNOTRemoved'+ ASType +'_Counts_ThresholdByPSI_' + str(PSI_Threshold))
    if FilterForeQTLS:
        Out = Out[Out['Hugo_Symbol'] == 0] # remove AS events that overlap with mutations
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE') # add cancer type to compare exp levels
    Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
    Out = Out.merge(Samples, left_on='Barcode', right_on='Barcode', how='left')
    Out['NMD'] = (Out['Group'] == 'Spliced_In') & (Out['DepletedTranscript'] == True)
    Out.groupby(['Barcode','MutLoad','NMD']).size().unstack().reset_index().replace(np.nan,0).to_csv(
        OutDir + 'TCGA_'+ ASType +'_Counts_ThresholdByPSI_' + str(PSI_Threshold) + 'FiltereQTLS_' + str(FilterForeQTLS))
    print('Done!')
    return(Out)



def GetDeltaPSIForEachGene():
    '''
    Performs a t-test for each intron retention event in TCGA between low and high mutational load tumors,
    removes intron retention events with missing calls and returns a table of average PSI changes for each gene,
    including if the change is significant.
    '''
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
    Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
    LowAndHigh = Samples[Samples['MutLoad'] < 10]['Barcode'].tolist() + Samples[Samples['MutLoad'] > 1000]['Barcode'].tolist()
    AS = GetAlternativeSplicingData(AS='RI', BarcodesOfInt=LowAndHigh)
    # Remove AS events if more than 25% of samples have a missing call
    PercentMissing = AS.isnull().sum(axis=1).reset_index()[0]/len(AS.columns)
    AS = AS.reset_index().assign(pct_missing = PercentMissing)
    AS = AS[AS['pct_missing'] < 0.25]
    # Get significantly different genes
    df = AS.apply(DoTTestForAS, axis=1).assign(GeneName=AS['symbol'])
    return(df)


def DoTTestForAS(row):
    '''
    For each AS event, performs a two-sided t-test for the null hypothesis that the two groups 
    (low and high mutational load tumors) have the same expected average PSI values. 
    Outputs a Pandas Series of p-values for how different the two groups are and the average 
    change in PSI values between the two groups for each gene.
    Parameters 
    @row = a series of PSI values for each patient from @GetDeltaPSIForEachGene 
    '''
    Samples = pd.read_csv(os.getcwd() + '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE')
    Samples['Barcode'] = Samples['Barcode'].str[0:12].str.replace('-','_')
    Low = row[list(set(row.index) & set(Samples[Samples['MutLoad'] < 10]['Barcode'].tolist()))].dropna()
    High = row[list(set(row.index) & set(Samples[Samples['MutLoad'] > 1000]['Barcode'].tolist()))].dropna()
    return( pd.Series({'pVal': stats.ttest_ind(Low, High)[1], 
						'Low_mean': Low.mean(), 
						'High_mean': High.mean(), 
						'LowToHighDelta': Low.mean() - High.mean()}))


def GetNumberGenesFilteredDueToPotentialEQTLs():
    '''
    Uses output from @GetExpLevelsForASEvents() to examine how many AS events are removed in
    both count categories (i.e., events splcied in/events not spliced in).
    '''
    ASDir='/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    Filt = pd.read_csv(ASDir +'TCGA_AS_EventsRemovedRI_Counts_ThresholdByPSI_0.8')
    Filt['Filtered'] = 1
    NotFilt= pd.read_csv(ASDir + 'TCGA_AS_EventsNOTRemovedRI_Counts_ThresholdByPSI_0.8')
    NotFilt['NotFiltered'] = 1
    FilteredStats = Filt.groupby(['Group'])['Filtered'].size().reset_index().merge(
        NotFilt.groupby(['Group'])['NotFiltered'].size().reset_index(),
        left_on=['Group'], right_on=['Group'])
    return(FilteredStats)



    
