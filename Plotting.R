PlotDir='/labs/ccurtis2/tilk/scripts/protein/Figures/'
library(scales)
library(cowplot)
library(dplyr)

AddPercentileRank = function(df) {
    df$PvaluePercentileRank = trunc(rank(df$pVal))/length(df$pVal)
    df$EstimatePercentileRankAll = trunc(rank(df$Estimate))/length(df$Estimate)
    out = rbind(
        data.frame(EstimatePercentileRank = trunc(rank(df[df['Estimate'] > 0,]$Estimate))/length(df[df['Estimate'] > 0,]$Estimate), 
        Sign = 'PositivePercentileRank', df[df['Estimate'] > 0,]), 
        data.frame(EstimatePercentileRank = trunc(rank(df[df['Estimate'] < 0,]$Estimate))/length(df[df['Estimate'] < 0,]$Estimate), 
        Sign = 'NegativePercentileRank', df[df['Estimate'] < 0,])
    )
    return(out)
}




PlotGlobalGSETCGA = function(df, Identifier) {
    # Match CORUM terms to broader categories
    df$grouping = ''
    df[grep('ribo',df[,1]),]$grouping = 'Translation' 
    df[grep('roteasome|PA',df[,1]),]$grouping = 'Proteasome' 
    df[grep('CCT',df[,1]),]$grouping = 'Chaperones'
    df[grep('synthesome|BRCA|FA|BLM|MCM|RC|BRAFT|DNA-PK',df[,1]),]$grouping = 'DNA Replication/Repair'
    df[grep('snRNP|plice|SMN|Rnase|Sm|Exosome',df[,1]),]$grouping = 'Splicing'
    df[grep('RC|CEN|Nup',df[,1]),]$grouping = 'Cell Cycle'
    # Add colors
    colors = data.frame(
        grouping = c('Transcription','Splicing','Chaperones','Translation','Cell Cycle', 'Chromosome Segregation','Proteasome', 'DNA Replication/Repair'),
        colors = as.character(c('#82CEF5','#3366C7','#A3E189','#3EA612','#F97FA2','#F73036','#F7B530','#FB791E'))
    )
    df = merge(df, colors, by='grouping')
    # Rank terms by significance for plotting 
    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    rankingDF = df %>% group_by(grouping) %>% summarise(max = max(negative_log10_of_adjusted_p_value))
    rankingDF = rankingDF[order(rankingDF$max),]
    rankingDF$rank = nrow(rankingDF):1
    df = merge(df, rankingDF, by='grouping')
    df = df[order(df$rank, -df$negative_log10_of_adjusted_p_value), ]
    df$ID = 1:nrow(df)
    df$ID2 = factor(df$ID, levels=df$ID)
    df$grouping = factor(df$grouping, levels=unique(df$grouping))
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
    PlotOut = ggplot(df, aes(x=reorder(term_name,-ID), y=negative_log10_of_adjusted_p_value, fill=colors)) +
        geom_col() +
        coord_flip() +
        facet_grid(rows=vars(grouping), scales='free_y', space='free', switch='y') +
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        #scale_fill_brewer(palette="Paired") +
        scale_fill_identity(guide = "legend", 
                            labels = unique(df$grouping),
                            breaks = unique(df$colors))  +
        # scale_fill_manual(values = df$colors, labels=df$grouping) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank())+
        theme(strip.text.y = element_blank())
    if (Identifier == 'TCGA_Regression') {
        ggsave(paste0(PlotDir, 'MixedEffectRregressionGeneSetEnrichmentGlobalTCGA.pdf' ), width=8, height=5, units='in') 
    } else if (Identifier == 'TCGA_DGE') {
         ggsave(paste0(PlotDir, 'DGEGeneSetEnrichmentGlobalTCGA.pdf' ), width=8, height=5, units='in') 
    } else if (Identifier == 'CCLE_Regression') {
        ggsave(paste0(PlotDir, 'RegressionGeneSetEnrichmentGlobalCCLE.pdf' ), width=8, height=5, units='in') 
    }
}




PlotRegCoefPerGroup = function(df, Dataset) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
    Estimate = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(Estimate)), color=subgroup)) +  geom_boxplot() +
                 theme_minimal() + labs(x='Effect Size (Beta Coefficient)', y='') + theme(legend.position='none') +
                 geom_vline(xintercept=0, linetype='dashed', col = 'black') 
    Rank = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(NegLog10Pval)), color=subgroup)) +  geom_boxplot() +
                theme_minimal() + labs(x='Negative Log10 of Adjusted P-Value', y='') + theme(legend.position='none') +
                theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
                geom_vline(xintercept= -log10(0.05), linetype='dashed', col = 'black')
    Combined = plot_grid(Estimate, Rank , ncol = 2) 
    if (Dataset == 'TCGA') {
        Title = ggdraw() + draw_label("TCGA", fontface='bold')
    } else if (Dataset == 'CCLE') {
        Title = ggdraw() + draw_label("CCLE", fontface='bold')
    } else if (Dataset == 'Protein_CCLE') {
        Title = ggdraw() + draw_label("CCLE", fontface='bold')
    }
    plot_grid(Title, Combined, ncol=1, rel_heights=c(0.1, 1)) 
    ggsave(paste0(PlotDir, 'RegCoefPerGroups', Dataset, '.pdf' ), width=8, height=5, units='in')
}



PlotStabilityPerTumor = function(df) {
    print(head(df))
    df$LogMutLoad = log10(as.numeric(as.character(df$MutLoad)))
    df$FractionUnstable = df$Unstable/(df$Stable + df$Unstable)
    ggplot(data=df, aes(x = LogMutLoad, y=FractionUnstable, color=Driver, fill=Driver)) +  geom_point() +
    facet_wrap(~Driver) + theme_minimal()
    ggsave(paste0(PlotDir, 'Stability_TCGA.pdf' ), width=8, height=8, units='in')
}

PlotDeltaGPerTumor = function(df) {
    print(head(df))
    df$LogMutLoad = log10(as.numeric(as.character(df$MutLoad)))
    ggplot(data=df, aes(x = LogMutLoad, y=(ddG), color=Driver, fill=Driver)) +  geom_point() +
    scale_y_continuous(limits= c(0.25,55), breaks = c(0.25,1,2,8,20,50), trans = log_trans()) +
    facet_wrap(~Driver) + theme_minimal()
    ggsave(paste0(PlotDir, 'Stability_TCGA.pdf' ), width=12, height=8, units='in')
}