PlotDir='/labs/ccurtis2/tilk/scripts/protein/Figures/'
library(scales)
library(cowplot)
library(dplyr)
library(ggpubr)
library(boot)


BootAS = function(df, NumBoot=100) {
    out = data.frame()
    for (i in 1:NumBoot) {
        rep = data.frame(df %>% group_by(Bin) %>% sample_frac(1, replace=TRUE) %>% 
                    summarize(retained = sum(True), not_retained = sum(False)) %>% 
                    group_by(Bin) %>% 
                    summarize(fraction = retained/(retained + not_retained)))
        rep$num = i
        out = rbind(out, rep)
    }
    return(out)
}

VisualizeAllAS = function(Threshold='0.8') {
    Dir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    AS_Types = c('RI','ES','AD','AP','AT','AA','ME')
    Combined = data.frame()
    df = data.frame()
    for (AS in AS_Types) {
        tmp = read.table(paste0(Dir, 'TCGA_',AS,'_Counts_ThresholdByPSI_', Threshold), sep=',', header=TRUE)
        tmp$AS = AS
        df = rbind(tmp, df)
    }
    df$Bin = cut(df$MutLoad, breaks=c(0,10,100,1000,50000), labels=c('0-10','10-100','100-1000','1000->10,000'))
    df = na.omit(df)
    df$FractionNMD= df$True/(df$False + df$True)
    booted= BootAS(df)
    booted_CI = data.frame(booted %>% group_by(Bin, AS) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    PlotOut = ggplot(booted_CI, aes(x=as.character(Bin), y=mean, group=AS, color=AS)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2) +
        labs(y='Under-expressed Transcripts With AS Event \nAll Transcripts With AS Event',
         x= 'Number of Protein Coding Mutations')

    ggsave(paste0(PlotDir, 'AS_RI_AllASType_TCGA.pdf' ), width=8, height=5, units='in') 


}






VisualizeAllASThresholds = function(AS_Type='RI') {
    Dir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    Thresholds = c('0.5','0.6','0.7','0.8','0.9')
    df = data.frame()
    for (T in Thresholds) {
        df = read.table(paste0(Dir, 'TCGA_',AS_Type,'_Counts_ThresholdByPSI_', T), sep=',', header=TRUE)
        df$Thresholds = T
        df = rbind(Combined, df)
    }
    df$Bin = cut(df$MutLoad, breaks=c(0,10,100,1000,50000), labels=c('0-10','10-100','100-1000','1000->10,000'))
    df = na.omit(df)
    df$FractionNMD= df$True/(df$False + df$True)
    booted= BootAS(df)
    booted_CI = data.frame(booted %>% group_by(Bin, Thresholds) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    PlotOut = ggplot(booted_CI, aes(x=as.character(Bin), y=mean, group=Thresholds, color=Thresholds)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2) +
        labs(y='Under-expressed Transcripts With Intron Retention \nAll Transcripts With Intron Retention',
         x= 'Number of Protein Coding Mutations')
    ggsave(paste0(PlotDir, 'AS_RI_AllThresholds_TCGA.pdf' ), width=8, height=5, units='in') 
}



VisualizeAS = function() {
    Dir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    df = read.table('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_RI_Counts_ThresholdByPSI_0.8' ,sep=',', header=T)
    df$Bin = cut(df$MutLoad, breaks=c(0,10,100,1000,50000), labels=c('0-10','10-100','100-1000','1000->10,000'))
    df = na.omit(df)
    df$FractionNMD= df$True/(df$False + df$True)
    foo = df %>%  group_by(Bin) %>% sample_frac(1, replace=TRUE) %>% 
                    summarize(retained = sum(True), not_retained = sum(False)) %>% 
                    group_by(Bin) %>% 
                    summarize(fraction = retained/(retained + not_retained))
    booted= BootAS(df)
    booted_CI = data.frame(booted %>% group_by(Bin) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    PlotOut = ggplot(booted_CI, aes(x=as.character(Bin), y=mean, group=1)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2, position=position_dodge(.9)) +
        labs(y='Under-expressed Transcripts With Intron Retention \nAll Transcripts With Intron Retention', x= 'Number of Protein Coding Mutations')
    ggsave(paste0(PlotDir, 'AS_RI_PerTumor_PSI_TCGA.pdf' ), width=5, height=4, units='in') 

}


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

PlotshRNA = function(df) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df$Group = gsub('mic Rib', 'mic\nRib', df$Group); df$Group = gsub('ial Rib', 'ial\nRib', df$Group); df$Group = gsub('ial Ch', 'ial\nCh', df$Group)
    df$Group = gsub('ory Pa', 'ory\nPa', df$Group)
    ggplot(df, aes(y=Dummy, x = Group ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","grey","blue"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5, space='free_x') +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='', fill='Beta\nCoefficient') + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=18))+
        ggtitle('Achilles (shRNA Screen)')
    ggsave(paste0(PlotDir, 'CCLE_shRNA_Regression.pdf' ), width=12, height=3.7, units='in')
}


PlotDrug = function(df) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df = subset(df, (df$subgroup != 'Protein Synthesis Inhibitor') & (df$subgroup != 'RNA Synthesis Inhibitor'))
    df$subgroup2 = gsub(' Inhibitor', '\nInhibitor', df$subgroup)
    df$subgroup2 = gsub('cific Pro', 'cific\nPro', df$subgroup)
    ggplot(df, aes(y=Dummy, x = name ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","white","grey"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5) +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup2), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='Drug Name', fill='Beta\nCoefficient') + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=14)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ggtitle('PRISM (Drug Screen)') 
    ggsave(paste0(PlotDir, 'CCLE_DrugRegression.pdf' ), width=12, height=2.75, units='in') 
}

PlotGlobalDownAndUpregulation= function(df) {
    all = subset(df, (df$value == 'Drivers') | (df$value == 'Passengers'))
    all$value = 'All'
    all = data.frame(all %>% group_by(value,Bin) %>% summarise(UnderExpressed = sum(UnderExpressed),NotUnderExpressed = sum(NotUnderExpressed),
             OverExpressed= sum(OverExpressed), NotOverExpressed= sum(NotOverExpressed)))
    df = rbind(all, df)
    #df$Bin2 = factor(df$Bin, levels= c('0 - 50', '50 - 100', '100 - 500', '500-1,000', '1,000-5,000', '5,000 - >10,000'))
    df$Bin2 = factor(df$Bin, levels= c('0 - 100', '100 - 1,000', '1,0000 - >10,000' ))
    df = subset(df, df$value %in% c('All','Drivers','Passengers','Essential','Housekeeping'))
    df$FractionOver = df$OverExpressed / (df$OverExpressed + df$NotOverExpressed)
    df$FractionUnder = df$UnderExpressed / (df$UnderExpressed + df$NotUnderExpressed)
    ggplot(df, aes(x=Bin2,y=FractionOver, color=value, group=value))+ geom_point() + geom_line() + 
        labs(y='Number Of Under-Expressed Genes\nTotal Number of Genes',x='Total Number of Protein Coding Mutations') +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        #facet_wrap(Essential~Drivers)
        ggsave(paste0(PlotDir, 'TCGA_UnderGlobalExpDriver.pdf' ), width=6, height=5, units='in') 

}

PlotGlobalKeggTCGA = function(df) {
    df = subset(df, df$source == 'KEGG')
    df = df[!grepl("disease", df$term_name),]
    df = df[!df$term_name %in% c('Amyotrophic lateral sclerosis','Epstein-Barr virus infection','Human T-cell leukemia virus 1 infection',
                            'Viral myocarditis','Influenza A','Spinocerebellar ataxia','Legionellosis','Biosynthesis of cofactors','Carbon metabolism'),]
    df$grouping = ''
    df[grep('Ribo|Pro|amino|Phagosome',df[,1]),]$grouping = 'Protein' 
    df[grep('RNA|Splic|mRNA',df[,1]),]$grouping = 'Transcription' 
    df[grep('repair|recombination|Cell|replication|sugar|backbone',df[,1]),]$grouping = 'DNA Replication/Repair'
    df[grep('Apop|senesc',df[,1]),]$grouping = 'Other'
    df[grep('Antigen|killer',df[,1]),]$grouping = 'Immune'
    colors = data.frame(
        grouping = c('Transcription','Protein','DNA Replication/Repair', 'Immune','Other' ),
        colors = as.character(c('#82CEF5','#3366C7','#3EA612','#F97FA2','#F73036'))
    )
    df = merge(df, colors, by='grouping')
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
        theme(strip.text.y = element_blank(), legend.position='bottom') +
        guides(fill=guide_legend(nrow=2,byrow=TRUE))
    ggsave(paste0(PlotDir, 'RegressionGeneSetEnrichmentGlobalKEGGTCGA.pdf' ), width=7, height=5, units='in') 

}

PlotGlobalGSETCGA = function(df, Identifier) {
    df= subset(df, df$source == 'CORUM')
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
        theme(strip.text.y = element_blank(), legend.position='bottom') + guides(fill=guide_legend(nrow=2,byrow=TRUE))
    if (Identifier == 'TCGA_Regression') {
        ggsave(paste0(PlotDir, 'MixedEffectRregressionGeneSetEnrichmentGlobalTCGA.pdf' ), width=7, height=5, units='in') 
    } else if (Identifier == 'TCGA_DGE') {
         ggsave(paste0(PlotDir, 'DGEGeneSetEnrichmentGlobalTCGA.pdf' ), width=7, height=5, units='in') 
    } else if (Identifier == 'CCLE_Regression') {
        ggsave(paste0(PlotDir, 'RegressionGeneSetEnrichmentGlobalCCLE.pdf' ), width=7, height=5, units='in') 
    }
}



PlotRefCoefAllGenes = function(df) {
    df$grouping = ''
    #df[grep('ribo',df[,2]),]$grouping = 'Translation' 
    df[grep('roteasome|PA',df[,2]),]$grouping = 'Proteasome' 
    df[grep('CCT',df[,2]),]$grouping = 'Chaperones'
    df[grep('synthesome|BRCA|FA|BLM|MCM|RC|BRAFT|DNA-PK',df[,2]),]$grouping = 'DNA Replication/Repair'
    df[grep('snRNP|plice|SMN|Rnase|Sm|Exosome',df[,2]),]$grouping = 'Splicing'
    df[grep('RC|CEN|Nup',df[,2]),]$grouping = 'Cell Cycle'
    # Add colors
    colors = data.frame(
        grouping = c('Transcription','Splicing','Chaperones','Translation','Cell Cycle', 'Chromosome Segregation','Proteasome', 'DNA Replication/Repair'),
        colors = as.character(c('#82CEF5','#3366C7','#A3E189','#3EA612','#F97FA2','#F73036','#F7B530','#FB791E'))
    )
    df = merge(df, colors, by='grouping')
    df$NegLog10PVal = as.numeric(as.character(df$NegLog10PVal))
    df$pval_rank = rank(df$NegLog10PVal, ties.method='min')
    ggplot(df, aes(x=NegLog10PVal, y=NegLog10PVal, fill=colors)) + geom_col() + theme_minimal() + 
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),legend.position = "none") +
        facet_wrap(~grouping, scales='free_x')
     ggsave(paste0(PlotDir, 'IndividualGeneFromGSEWithGlobalExpressionTCGA.pdf' ), width=8, height=5, units='in')

}

PlotRegCoefPerGroup = function(df) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$Dataset = gsub('TCGA', 'TCGA (Human Tumors)', df$Dataset); df$Dataset = gsub('CCLE', 'CCLE (Cancer Cell Lines)', df$Dataset)
    df$Dataset2 = factor(df$Dataset, levels=c('TCGA (Human Tumors)','CCLE (Cancer Cell Lines)'))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
    Estimate = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(Estimate)), color=subgroup)) +  geom_boxplot() +
                 theme_minimal() + labs(x='Effect Size (Beta Coefficient)', y='') + theme(legend.position='none') +
                 geom_vline(xintercept=0, linetype='dashed', col = 'black') + facet_wrap(~Dataset2, scales='free_x') +
                 theme(strip.text = element_text(face="bold", size=12))

    Rank = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(NegLog10Pval)), color=subgroup)) +  geom_boxplot() +
                theme_minimal() + labs(x='Negative Log10 of Adjusted P-Value', y='') + 
                theme(legend.position='bottom', legend.title=element_blank()) +
                #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
                geom_vline(xintercept= -log10(0.05), linetype='dashed', col = 'black') + facet_wrap(~Dataset2, scales='free_x') +
                theme(strip.text = element_text(face="bold", size=12))
    Combined = plot_grid(Estimate, Rank, rel_heights=c(0.85, 1), ncol = 1) 
    #plot_grid(Title, Combined, ncol=1, rel_heights=c(0.1, 1)) 
    ggsave(paste0(PlotDir, 'RegCoefPerGroups_TCGAandCCCLE.pdf' ), width=6, height=7, units='in')
}

PlotDeltaPSI = function(df) {
    print(df)
    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    df$Group = gsub('PosPSI','Positive Delta PSI\n(Low to High)', df$Group)
    df$Group = gsub('NegPSI','Negative Delta PSI\n(Low to High)', df$Group)
    rankingDF = df %>% group_by(Group) %>% summarise(max = max(negative_log10_of_adjusted_p_value))
    rankingDF = rankingDF[order(rankingDF$max),]
    rankingDF$rank = nrow(rankingDF):1
    df = merge(df, rankingDF, by='Group')
    df = df[order(df$rank, -df$negative_log10_of_adjusted_p_value), ]
    df$ID = 1:nrow(df)
    df$ID2 = factor(df$ID, levels=df$ID)
    df$grouping = factor(df$Group, levels=unique(df$Group))
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
    Out = ggplot(df, aes(x=reorder(term_name,-ID), y=negative_log10_of_adjusted_p_value, fill=source)) +
        geom_col() +
        coord_flip() +
        facet_grid(rows=vars(Group), scales='free_y', space='free', switch='y') +
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        scale_fill_manual(values = c('red','orange')) +
        theme_minimal() +
        theme( legend.title=element_blank())
    ggsave(paste0(PlotDir, 'AS_Delta_PSI_IntonRetention_TCGA.pdf' ), width=7, height=4.5, units='in')
}

PlotProtein = function(df) {
    df$Value = as.numeric(as.character(df$Value))
    df = df[df$Group %in% c('Cytoplasmic Ribosomes','Mitochondrial Ribosomes'),]
    ggviolin(df, x = "MutGroup", y = "Value", xlab = "", ylab = "Normalized Protein Expression",
        add = "boxplot", add.params = list(fill = "white"),
        fill = "Group", palette = "jco", facet.by = "Group", ylim = c(-2.5,2.5)) + 
        stat_compare_means(method = "t.test", label =  "p.signif", label.x = 1.5)
    # ggplot(data=df, aes(x = MutGroup, y=Value, color=Group)) +  geom_boxplot() + 
    # theme_minimal() + labs(y='Normalized Protein Expression',x='') + facet_wrsap(~Group)
    ggsave(paste0(PlotDir, 'Protein_Expression_CCLE.pdf' ), width=5, height=5, units='in')

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


PlotOverAndUnderForProteasome = function() {

    library(reshape2)
    df = read.table('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/NumOverAndOverForProteasomeByGene', sep=',', header=T, stringsAsFactors=FALSE)
    df = df[df$MutBin %in% c('high','low'),]
    df$X = NULL 
    
    plot_df = melt(df[c('MutBin','group','FractionOfGenesInCategory')], id.vars = c("group",'MutBin'))
    #plot_df = subset(plot_df, plot_df$variable == 'DepletedTranscript')
    PlotOut = ggplot(df, aes(x=MutBin, y=FractionOfGenesInCategory, color=group)) + geom_boxplot() + facet_wrap(~type)
        #facet_wrap(~MutLoadGroup, scale='free_y')
    ggsave(paste0(PlotDir, 'ProteasomeExpressionForClassification.pdf' ), width=8, height=5, units='in')

}

PlotJacknifedExpressionAcrossGroups = function(df, dataset) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
    Estimate = ggplot(data=df, aes(x = CancerTypeRemoved, y=as.numeric(as.character(Estimate)), color=SortedLevel)) +  geom_boxplot() +
                facet_wrap(~SortedLevel, ncol=5) +
                 theme_minimal() + labs(y='Effect Size (Beta Coefficient)', x='') +  theme(legend.position='none') +
                 geom_hline(yintercept=0, linetype='dashed', col = 'black') + #facet_wrap(~SortedLevel, scales='free_x', ncol=12) +
                 theme(strip.text = element_text(face="bold", size=8)) + ggtitle(dataset) +
                 theme(plot.title = element_text(hjust = 0.5)) + # Center title
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))
    ggsave(paste0(PlotDir, 'RegCoefPerGroupsJacknifed_', dataset,'.pdf' ), width=10, height=6, units='in')

}

PlotJacknifedDrugsAcrossGroups = function(df) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df = subset(df, (df$subgroup != 'Protein Synthesis Inhibitor') & (df$subgroup != 'RNA Synthesis Inhibitor'))
    df$subgroup2 = gsub(' Inhibitor', '\nInhibitor', df$subgroup)
    df$subgroup2 = gsub('cific Pro', 'cific\nPro', df$subgroup)
    ggplot(df, aes(y=CancerTypeRemoved, x = name ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","white","grey"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5) +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup2), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='Drug Name', fill='Beta\nCoefficient') + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=12)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ggtitle('PRISM (Drug Screen)') 
    ggsave(paste0(PlotDir, 'CCLEDrugJacknifed.pdf' ), width=12, height=6, units='in')
}


PlotJacknifedshRNAAcrossGroups = function(df) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df$Group = gsub('mic Rib', 'mic\nRib', df$Group); df$Group = gsub('ial Rib', 'ial\nRib', df$Group); df$Group = gsub('ial Ch', 'ial\nCh', df$Group)
    df$Group = gsub('ory Pa', 'ory\nPa', df$Group)
    ggplot(df, aes(y=CancerTypeRemoved, x = Group ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","grey","blue"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5, space='free_x') +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='', fill='Beta\nCoefficient') + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=18))+
        ggtitle('Achilles (shRNA Screen)')
    ggsave(paste0(PlotDir, 'CCLEshRNAJacknifed.pdf' ), width=12, height=8, units='in')
}


PlotRegressionForAllDrugs = function(df) {
    df$log10pval= log10(as.numeric(as.character(df$pVal)))
    df$SigGroups = ifelse(df$pVal < 0.05, 'Significant', 'NotSignificant')
    df$EstimateGroups = ifelse(df$Estimate > 0, 'Increase In Viability With Load','Decrease In Viability With Load')
    #print(head(df))
    PlotOfCounts = ggplot(df, aes(x = EstimateGroups, color=SigGroups, fill=SigGroups)) + geom_bar() +
        theme_minimal() + labs(y='Number of Drugs', x='') + 
        theme(legend.position='bottom', legend.title=element_blank()) 
    df = subset(df, (df$SigGroups == 'Significant'))
    CountsPerDrug = data.frame(table(as.character(df$subgroup)))
    CountsPerDrug$NewGroup = ifelse(CountsPerDrug$Freq < 2, 'Other', as.character(CountsPerDrug$Var1))
    df= merge(df, CountsPerDrug, by.x='subgroup', by.y='Var1') 
    df = subset(df, (df$Estimate < 0))
    PlottingGroup = rbind(
        data.frame(PlottingGroup = 'Other', NewGroup = c('Other')),
        data.frame(PlottingGroup = 'Growth Factor Inhibitors', NewGroup = c("acetylcholine receptor antagonist","adenosine receptor antagonist",
                                        "adrenergic receptor agonist","adrenergic receptor antagonist",
                                        'benzodiazepine receptor agonist','androgen receptor antagonist','dipeptidyl peptidase inhibitor',
                                        'dopamine receptor agonist','dopamine receptor antagonist','dopamine receptor antagonist, serotonin receptor antagonist',
                                        'EGFR inhibitor',"FGFR inhibitor, VEGFR inhibitor","glutamate receptor antagonist","PDGFR tyrosine kinase receptor inhibitor",
                                        "progesterone receptor agonist","retinoid receptor agonist",'serotonin receptor agonist','serotonin receptor antagonist','VEGFR inhibitor')),
        data.frame(PlottingGroup = 'Nuclear Transport Inhibitor', NewGroup = c('exportin antagonist')),
        data.frame(PlottingGroup = 'Pro-Apoptosis' , NewGroup = c("AKT inhibitor","ALK tyrosine kinase receptor inhibitor","Aurora kinase inhibitor",'CDK inhibitor',
                                        'JNK inhibitor','mTOR inhibitor, PI3K inhibitor','NFkB pathway inhibitor','other antibiotic','phosphodiesterase inhibitor',
                                        'PI3K inhibitor','PKC inhibitor','PLK inhibitor','protein tyrosine kinase inhibitor','RAF inhibitor','tyrosine kinase inhibitor','XIAP inhibitor')),
        data.frame(PlottingGroup = 'ATPase Inhibitor' , NewGroup = c("ATPase inhibitor",'HDAC inhibitor')),
        data.frame(PlottingGroup = 'Inflammatory/Immune' , NewGroup = c('analgesic agent','CC chemokine receptor antagonist','cyclooxygenase inhibitor',
                                        'glucocorticoid receptor agonist','histamine receptor antagonist','nitric oxide synthase inhibitor','local anesthetic',
                                        'p38 MAPK inhibitor','SYK inhibitor')),
        data.frame(PlottingGroup = 'Protein Synthesis Inhibitor' , NewGroup = c('antimalarial agent','bacterial 50S ribosomal subunit inhibitor','protein synthesis inhibitor')),
        data.frame(PlottingGroup = 'DNA Replication Inhibitor' , NewGroup = c('antiprotozoal agent','ATR kinase inhibitor','bacterial DNA gyrase inhibitor','CHK inhibitor',
                                        'DNA inhibitor','DNA polymerase inhibitor','ribonucleotide reductase inhibitor','topoisomerase inhibitor')),
        data.frame(PlottingGroup = 'Microtubulin/Spindle Polymerization Inhibitor' , NewGroup = c('bacterial cell wall synthesis inhibitor','fungal lanosterol demethylase inhibitor',
                                        'kinesin-like spindle protein inhibitor','membrane integrity inhibitor','microtubule inhibitor',
                                        'microtubule stabilizing agent, tubulin polymerization inhibitor','tubulin polymerization inhibitor')),
        data.frame(PlottingGroup = 'Transcription Inhibitor' , NewGroup = c('bromodomain inhibitor','RNA polymerase inhibitor')),
        data.frame(PlottingGroup = 'Ion Channel Regulation', NewGroup =  c('calcium channel blocker','potassium channel blocker','potassium channel activator','sodium channel blocker')),
        data.frame(PlottingGroup = 'Cell Migration Inhibitor', NewGroup =  c("focal adhesion kinase inhibitor",'matrix metalloprotease inhibito')),
        data.frame(PlottingGroup = 'Protein Degradation Inhibitors', NewGroup =  c("HCV inhibitor",'HIV protease inhibitor',"proteasome inhibitor",'ubiquitin specific protease inhibitor')),
        data.frame(PlottingGroup = 'Chaperone Inhibitors', NewGroup =  c("HSP inhibitor")),
        data.frame(PlottingGroup = 'Oxidative Stress Alleviation Inhibitor', NewGroup =  c("nuclear factor erythroid derived, like (NRF2) activator"))
    )
    df = merge(df, PlottingGroup, by.x='NewGroup' ,by.y='NewGroup')
    df$EstimateRank = rank(df$Estimate)
    PlotFactors = data.frame(df %>% group_by(PlottingGroup) %>% summarize(median(Estimate)))
    PlotFactors = PlotFactors[order(PlotFactors$median.Estimate.),]
    df$PlottingGroups2 = factor(df$PlottingGroup, levels=as.character(PlotFactors$PlottingGroup))
    Foo = ggplot(data=df, aes(x = Estimate, y=PlottingGroups2, color=PlottingGroups2)) +  geom_boxplot() +
        theme_minimal() + labs(x='Effect Size (Beta Coefficient)', y='') +
        #theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        geom_hline(yintercept= 0, linetype='dashed', col = 'black') +theme(legend.position='none') +
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
    Combined = plot_grid(PlotOfCounts, Foo, rel_heights=c(0.75,1), ncol = 1) 
    ggsave(paste0(PlotDir, 'RegCoefAllDrugs_CCLE.pdf' ), width=6, height=6, units='in')
}
