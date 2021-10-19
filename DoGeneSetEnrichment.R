library('gprofiler2')


ReadDGEResults = function(Category) {
    Dir="/labs/ccurtis2/tilk/scripts/protein/Data/DGE/Output/"
    if (Category == 'PanCancerDGE') {
        return(read.table(paste0(Dir, 'TCGA_PanCancer_DGE.txt'), header=TRUE))
    } else if (Category == 'PanCancerMatchedDGE') {
        return(read.table(paste0(Dir, 'TCGA_PanCancerMatched_DGE.txt'), header=TRUE))
    }
}

ReadRegressionResults = function(Category) {
    Dir="/labs/ccurtis2/tilk/scripts/protein/Data/Regression/"
    if (Category == 'TCGA') {
        df = read.table(paste0(Dir, 'tcga_expression_whole_genome.txt'), sep=',', header=TRUE)
    } else if (Category == 'CCLE') {
        df = read.table(paste0(Dir, 'ExpressionOLSRegressionEstimatesCCLE'), sep=',', header=TRUE)
        df$adj.pval = p.adjust(df$Pr...t.., method= 'fdr')
    } else if (Category == 'RNAiCCLE') {
        df = read.table(paste0(Dir, 'RNAiOLSRegressionEstimatesCCLE'), sep=',', header=TRUE)
        df$adj.pval = p.adjust(df$Pr...t.., method= 'fdr')
    }
    df$Coefficient = gsub('[[:digit:]]+', '', df$X)
    df = df[df['Coefficient'] == 'LogScore',]
    return(df)
}

GetGeneSets = function(Category) {
    NumberGenesToSubset = 5111
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
    return(gost(GeneSet)$result)
}




foo = DoGeneSetEnrichment(GetGeneSets('PanCancerMatchedDGE'))

foo = DoGeneSetEnrichment(GetGeneSets('PanCancerDGE'))

foo = DoGeneSetEnrichment(GetGeneSets('ExpressionTCGA'))

foo = DoGeneSetEnrichment(GetGeneSets('ExpressionCCLE'))


foo = DoGeneSetEnrichment(GetOverlapOfGeneSets)

foo[foo['source'] == 'CORUM',]$term_name

foo[foo['source'] == 'REAC',]$term_name