

ReadDGEResults = function(Category) {
    Dir="/labs/ccurtis2/tilk/scripts/protein/Data/DGE/Output/"
    if (Category == 'PanCancerDGE') {
        return(read.table(paste0(Dir, 'TCGA_PanCancer_DGE.txt'), header=TRUE))
    } else if (Category == 'PanCancerMatchedDGE') {
        return(read.table(paste0(Dir, 'TCGA_PanCancerMatched_DGE.txt'), header=TRUE))
    }
}

AddPercentileRank = function(df) {
    df$PvaluePercentileRank = trunc(rank(df$pVal))/length(df$pVal)
    out = rbind(
        data.frame(EstimatePercentileRank = trunc(rank(df[df['Estimate'] > 0,]$Estimate))/length(df[df['Estimate'] > 0,]$Estimate), 
        Sign = 'PositivePercentileRank', df[df['Estimate'] > 0,]), 
        data.frame(EstimatePercentileRank = trunc(rank(df[df['Estimate'] < 0,]$Estimate))/length(df[df['Estimate'] < 0,]$Estimate), 
        Sign = 'NegativePercentileRank', df[df['Estimate'] < 0,])
    )
    return(out)
}



ReadRegressionResults = function(Category) {
    Dir="/labs/ccurtis2/tilk/scripts/protein/Data/Regression/"
    if (Category == 'TCGA') {
       # df = read.table(paste0(Dir, 'ExpressionMixedEffectRegressionEstimatesKsKaTCGA'), sep=',', header=TRUE)
        df = read.table(paste0(Dir, 'ExpressionMixedEffectRegressionEstimatesKsKaTCGAPurity'), sep=',', header=TRUE)
    } else if (Category == 'CCLE') {
        #df = read.table(paste0(Dir, 'ExpressionOLSRegressionEstimatesCCLE'), sep=',', header=TRUE)
        df = read.table(paste0(Dir, 'ExpressionOLSRegressionEstimatesKsKaCCLE'), sep=',', header=TRUE)
        return(df)
    } else if (Category == 'RNAiCCLE') {
        df = read.table(paste0(Dir, 'RNAiOLSRegressionEstimatesKsKaCCLE'), sep=',', header=TRUE) 
        #df = read.table(paste0(Dir, 'StandardizedRNAiOLSRegressionEstimatesKsKaCCLE'), sep=',', header=TRUE)
    }
    # Calculate adjusted p-values for multiple hypothesis testing
    df$Coefficient = gsub('[[:digit:]]+', '', df$X)
    df = df[df['Coefficient'] == 'LogScore',]
    df$adj.pval = p.adjust(df$Pr...t.., method= 'fdr')
    return(df)
}

GetGeneSets = function(Category) {
    NumberGenesToSubset = 1000000
    if (Category == 'PanCancerDGE') {
        df = ReadDGEResults(Category)
        # Choose only significant genes that are upregulated in high TMB tumors
        GeneSet = subset(df, (df$padj < 0.05) & (df$log2FoldChange < 0))
        # This filter produces ~13k sig genes, filter them down to top # of genes
        GeneSet = row.names(head(GeneSet[order(GeneSet$log2FoldChange),], n=NumberGenesToSubset))
    } else if (Category == 'PanCancerMatchedDGE') {
         df = ReadDGEResults(Category)
         GeneSet = subset(df, (df$padj < 0.05) & (df$log2FoldChange < 0))
         GeneSet = row.names(head(GeneSet[order(GeneSet$log2FoldChange),], n=NumberGenesToSubset))
    } else if (Category == 'ExpressionTCGA') {
        df = ReadRegressionResults('TCGA')
        GeneSet = df[(df['Estimate'] > 0) & (df['adj.pval'] < 0.05),]$GeneName
    } else if (Category == 'ExpressionTCGALow') {
        df = ReadRegressionResults('TCGA')
        GeneSet = df[(df['Estimate'] < 0) & (df['adj.pval'] < 0.05),]$GeneName
    } else if (Category == 'ExpressionCCLE') {
        df = ReadRegressionResults('CCLE')
        GeneSet = df[(df['Estimate'] > 0) & (df['adj.pval'] < 0.05),]$GeneName
    } else if (Category == 'RNAiCCLE') {
        df = ReadRegressionResults('RNAiCCLE')
        GeneSet = df[(df['Estimate'] < 0) & (df['adj.pval'] < 0.05),]$GeneName
    }
    return(GeneSet)
}

GetOverlapOfGeneSets = function() {
    ccle = GetGeneSets('ExpressionCCLE')
    tcga = GetGeneSets('ExpressionTCGA')
    return(tcga[tcga %in% ccle])
}


DoGeneSetEnrichment = function(GeneSet) {
    library('gprofiler2')
    result = data.frame(gost(GeneSet)$result)
    return(result[c('source','term_name','p_value')])
}




# foo = DoGeneSetEnrichment(GetGeneSets('PanCancerMatchedDGE'))

# foo = DoGeneSetEnrichment(GetGeneSets('PanCancerDGE'))

# foo = DoGeneSetEnrichment(GetGeneSets('ExpressionTCGA'))
# foo = DoGeneSetEnrichment(GetGeneSets('ExpressionTCGA_Low'))


# foo = DoGeneSetEnrichment(GetGeneSets('RNAiCCLE'))

# foo = DoGeneSetEnrichment(GetGeneSets('ExpressionCCLE'))


# foo = DoGeneSetEnrichment(GetGeneSets('RNAiCCLE'))


# foo = DoGeneSetEnrichment(GetOverlapOfGeneSets())

# foo[foo['source'] == 'CORUM',]$term_name

# foo[foo['source'] == 'REAC',]$term_name